#include "fabm_driver.h"

module FEISTY_FABM

   use fabm_types                         ! FABM
   use fabm_particle                      ! FABM
   use fabm_builtin_depth_mapping         ! FABM
   use fabm_expressions                   ! FABM
   use feisty_fabm_vertical_distribution  ! FEISTY vertical distribution weights (feisty_layer_weights)

   use globals                            ! FEISTY
   use spectrum                           ! FEISTY
   use fish                               ! FEISTY
   use setup                              ! FEISTY
   
   implicit none

   private

   ! Fluxes from the fish to the BGC model: respiration, feces and carcasses, each as C, N and P
   integer, parameter :: n_flux = 9
   integer, parameter :: flux_resp_c = 1, flux_resp_n = 2, flux_resp_p = 3, &
                         flux_feces_c = 4, flux_feces_n = 5, flux_feces_p = 6, &
                         flux_carcass_c = 7, flux_carcass_n = 8, flux_carcass_p = 9
   character(len=13), parameter :: flux_names(n_flux) = [character(len=13) :: 'respiration_c', 'respiration_n', 'respiration_p', &
                                   'feces_c', 'feces_n', 'feces_p', 'carcasses_c', 'carcasses_n', 'carcasses_p']

   type, extends(type_depth_integrated_particle), public :: type_feisty_fabm
      ! State variables for fish (small pelagics, mesopelagics, large pelagics, mid-water predators, and demersal fish):
      type (type_bottom_state_variable_id),         allocatable :: id_fish(:)
      ! Uniform vertical distribution used for the depth-integrated zooplankton prey:
      type (type_vertical_distribution_id)                      :: dummy_w
      ! The benthos state variable:
      type (type_bottom_state_variable_id)                      :: id_benthos, id_det_bot_c, id_det_bot_n, id_det_bot_p, id_det_bot_s
      type (type_bottom_dependency_id),             allocatable :: id_ben_c(:), id_ben_n(:), id_ben_p(:)
      !type (type_bottom_state_variable_id)                      :: id_large_benthos
      type (type_state_variable_id)                             :: id_nut_c, id_nut_n, id_nut_p

      ! Dependency IDs for external FABM variables (e.g., small and large zooplankton, nutrients, detritus):
      type (type_bottom_dependency_id),    allocatable          :: id_smzoo_c(:), id_smzoo_n(:), id_smzoo_p(:) 
      type (type_bottom_dependency_id),    allocatable          :: id_lgzoo_c(:), id_lgzoo_n(:), id_lgzoo_p(:) 
      type (type_bottom_dependency_id)                          :: id_zoo_c, id_zoo_n, id_zoo_p, id_det_flux_c, id_det_flux_n, id_det_flux_p, id_det_flux_s
      ! Fluxes of each fish size class to the BGC model (fish, flux type), set in do_bottom and read in do_column:
      type (type_bottom_diagnostic_variable_id),    allocatable :: id_flux(:,:)
      type (type_horizontal_dependency_id),         allocatable :: id_flux_in(:,:)
      ! Pelagic pools receiving the vertically distributed fluxes, one per flux type:
      type (type_state_variable_id)                             :: id_flux_target(n_flux)
      logical                                                   :: flux_on(n_flux)
      ! Layer depth and thickness, for the vertical distribution of the fluxes:
      type (type_dependency_id)                                 :: id_z, id_h
      ! Coupling; pointer to which model contains the zooplankton state variable that we will integrater over the water column:
      type (type_model_id),    allocatable                      :: id_smzoo_int(:), id_benprey(:) 
      type (type_model_id),    allocatable                      :: id_lgzoo_int(:) 
      type (type_model_id)                                      :: id_zoo_int, id_det_target_c, id_det_target_n, id_det_target_p, id_det_target_s
      
      type (type_dependency_id)                                 :: id_temp
      type (type_horizontal_dependency_id)                      :: id_temp_vertmean_100m,id_temp_vertmean_500to1500m,id_bottom_depth

      type (type_bottom_diagnostic_variable_id)                 :: id_temp_vertmean_100m_diag, id_det_pool_c, id_det_pool_n, id_det_pool_p, id_det_pool_s
      type (type_bottom_diagnostic_variable_id), allocatable    :: id_fish_total_biomass(:)
      type (type_bottom_diagnostic_variable_id), allocatable    :: id_fish_catch(:), id_fish_F(:)
      ! Fishing mortality source per functional type (0: Fmax_i parameter, 1: external maximum F via dependency F_<fft>)
      integer,                                   allocatable    :: Fsource(:)
      type (type_horizontal_dependency_id),      allocatable    :: id_Fext(:)

      real (rk)                              :: b_ini
      real (rk) ,    allocatable             :: ben_availability(:), depth_array(:) ! Array of depth values for pre-computed profiles
      real (rk) ,    allocatable             :: zoo_availability(:), smzoo_availability(:), lgzoo_availability(:) ! availability of each zooplankton prey (n_zoo > 1): all, small, large
      real (rk) ,    allocatable             :: theta_matrix(:,:,:) ! Pre-computed theta matrix (predator,prey,n_depth)
      real (rk) ,    allocatable             :: vertical_distribution_matrix_centers(:,:,:) ! Pre-computed vertical distribution matrix (depth_center,nFGrid,n_depth) - contains center values
      
      real(rk)                               :: shelfdepth
      real(rk)                               :: max_depth
      real(rk)                               :: photic_depth, benk, ben_eff
      real(rk)                               :: det_c_scale, det_n_scale, det_p_scale, det_s_scale
      integer                                :: n_depth, n_zoo, n_s, n_l, n_ben=0
      character,    allocatable              :: zoo_size(:)
      character(len=100), allocatable        :: fft_short_names(:), size_numstr(:)
      character(len=4)                       :: respiration_elements, feces_elements, carcass_elements, bottom_detritus_elements
      real(rk)                               :: bgmort    ! Fish background mortality (per year)
      real(rk)                               :: gww_gc    ! gram wet weight to gram carbon conversion factor
      real(rk)                               :: qnc, qpc  ! Nitrogen:Carbon and Phosphorus:Carbon ratios
      !real(rk)                               :: mol_targetunit  ! mol to target unit conversion factor NIAZ: always the same when we couple through standard_variables, so we use a hardcoded conversion
      real(rk)                               :: gww_mmolc, gww_FISHn, gww_FISHp   ! gww to target C/N/P unit conversion factors

      logical                                :: internal_benthos !, nutrient_in_c, nutrient_in_n, nutrient_in_p           ! BGC has C/N/P elements for nutrient coupling
      logical                                :: resp_c, resp_n, resp_p, feces_c, feces_n, feces_p, carcass_c, carcass_n, carcass_p, bot_det_c, bot_det_n, bot_det_p, bot_det_s
      
   contains
      procedure :: initialize
      procedure :: do_bottom
      procedure :: do_column
      ! Reference model procedures here.
   end type type_feisty_fabm
   
   real(rk), parameter :: seconds_per_year = 365._rk*86400._rk

contains

   subroutine initialize(self, configunit)
      class (type_feisty_fabm), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
      
      real(rk)           :: szprod, lzprod, bprodin, dfbot, depth, Tp, Tb, bent_ini
      real(rk)           :: dfpho, Tm, etamature, visual, Fmax, etaF!, ssigma, tau, shelfdepth
      real(rk)           :: Fmax_i, etaF_i   ! per functional type fishing parameters
      integer            :: nStages
      integer            :: i,j,k, n_names
      character(len=100) :: i_str, j_str, size_number_str, fft_short_name
      character(len=200) :: flux_name

      class (type_vertical_depth_range), pointer :: depth_distribution

      !call self%type_depth_integrated_particle%initialize(configunit)
      call self%register_implemented_routines((/source_do_bottom, source_do_column/))
      
      ! Register model parameters  
      ! Fish physiological parameters. See input.nml in the original FEISTY R package and Fortran library
      call self%get_parameter(h, 'h', 'g^nn yr-1', 'Max. consumption coefficient', default=20._rk)  
      call self%get_parameter(nn, 'nn', '-', 'Metabolic exponent', default=-0.25_rk)  
      call self%get_parameter(gamma, 'gamma', 'm2 g^q yr-1', 'Coef. for clearance rate', default=70._rk)  
      call self%get_parameter(q, 'q', '-', 'Clearance rate exponent', default=-0.2_rk)  
      !call self%get_parameter(kk, 'kk', 'g^p yr-1', 'Metabolism coefficient', default=4._rk)  
      call self%get_parameter(p, 'p', '-', 'Metabolism exponent', default=-0.175_rk)  
      call self%get_parameter(epsAssim, 'epsAssim', '-', 'Assimilation efficiency', default=0.7_rk)
      call self%get_parameter(epsRepro, 'epsRepro', '-', 'Reproduction & recruitment efficiency', default=0.01_rk)  
      call self%get_parameter(self%bgmort, 'bgmort', 'yr-1', 'Fish background mortality', default=0.1_rk, scale_factor=1._rk/seconds_per_year) 
      
      call self%get_parameter(beta, 'beta', '-', 'Beta parameter for size-based predation preference', default=400._rk)  
      call self%get_parameter(sigma, 'sigma', '-', 'Sigma parameter for size-based predation preference', default=1.3_rk)  
      call self%get_parameter(mMin, 'mMin', 'g', 'Minimum fish mass (boundary of the grid)', default=0.001_rk)  
      !call self%get_parameter(mMedium, 'mMedium', 'g', 'Medium fish central mass for feeding preference calc', default=0.5_rk)   ! for vertical distribution 10
      !call self%get_parameter(mLarge, 'mLarge', 'g', 'Large fish central mass for feeding preference calc', default=250._rk)     ! for vertical distribution 5000
      
          
      !call self%get_parameter(lbenk, 'lbenk', 'g m-2', 'Large benthos carry capacity', default=0._rk)  !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(szoog, 'szoog', 'yr-1', 'Small zooplankton growth rate', default=1._rk)  !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(lzoog, 'lzoog', 'yr-1', 'Large zooplankton growth rate', default=1._rk)  !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(sbeng, 'sbeng', 'yr-1', 'Small benthos growth rate', default=1._rk)      !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(lbeng, 'lbeng', 'yr-1', 'Large benthos growth rate', default=0._rk)      !NIAZ: deleted (no zoop parameters should be in the fish code)
      ! predation preference coefficient
      !call self%get_parameter(thetaS, 'thetaS', '-', 'Medium fish preference for small zooplankton', default=0.25_rk)  
      !call self%get_parameter(thetaA, 'thetaA', '-', 'Large fish preference for medium forage fish', default=0.5_rk)  
      !call self%get_parameter(thetaD, 'thetaD', '-', 'Preference of large demersal on pelagic prey', default=0.75_rk)  
      ! get initial values
      !call self%get_parameter(smz_ini,'smz_ini', 'g m-2',  'initial small mesozooplankton biomass', default=1.e2_rk)  
      !call self%get_parameter(lgz_ini, 'lgz_ini', 'g m-2', 'initial large mesozooplankton biomass', default=1.e2_rk)  
      call self%get_parameter(bent_ini, 'bent_ini', 'g m-2', 'initial benthos biomass', default=5._rk)  
      !call self%get_parameter(lgbent_ini, 'lgbent_ini', 'g m-2', 'initial large benthos biomass', default=0._rk)        
      call self%get_parameter(self%b_ini, 'b_ini', 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  
      ! get input values
      !call self%get_parameter(szprod, 'szprod', 'g m-2 yr-1', 'small mesozooplankton production', default=100._rk)  
      !call self%get_parameter(lzprod, 'lzprod', 'g m-2 yr-1', 'large mesozooplankton production', default=100._rk)  
      !call self%get_parameter(bprodin, 'bprodin', 'g m-2 yr-1', 'benthos production', default=5._rk)  
      !call self%get_parameter(dfbot, 'dfbot', 'g m-2 yr-1', 'detrital flux reaching the bottom', default=-1._rk)  
      call self%get_parameter(depth, 'depth', 'm', 'water column depth', default=100._rk)  
      call self%get_parameter(Tp, 'Tp', 'Celsius', 'pelagic layer averaged temperature', default=10._rk)  
      call self%get_parameter(Tb, 'Tb', 'Celsius', 'bottom layer depth temperature (500m - up to 1500m)', default=8._rk) 
      
      !call self%get_parameter(dfpho, 'dfpho', 'g m-2 yr-1', 'detrital flux out of the photic zone', default=-1._rk)
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
      
      !call self%get_parameter(bET_val, 'bET_val', '-', 'whether turn on the effective temperature (integer 1 or 0)', default=1)
      
      call self%get_parameter(self%max_depth, 'max_depth', 'm', 'maximum depth of the model domain', default=5000._rk)
      call self%get_parameter(self%n_depth, 'n_depth', '-', 'number of depth levels for vertical distribution matrices', default=50)
      call self%get_parameter(self%gww_gc, 'gww_gc', '-', 'gram wet weight to gram carbon ratio', default=9._rk)
      call self%get_parameter(self%qnc, 'qnc', 'N/C', 'Nitrogen to Carbon ratio', default=16.0_rk/106.0_rk)
      call self%get_parameter(self%qpc, 'qpc', 'P/C', 'Phosphorus to Carbon ratio', default=1.0_rk/106.0_rk)
      !call self%get_parameter(self%mol_targetunit, 'mol_targetunit', '-', 'mol to target unit (in LTL model) ratio', default=1._rk/1000._rk)! default is mol m-3:mmol m-3 (used in LTL model)
      
      call self%get_parameter(self%feces_elements, 'feces_elements', '', 'elemental composition of feces e.g. "cnp"')
      self%feces_c = index(self%feces_elements, 'c') /= 0
      self%feces_n = index(self%feces_elements, 'n') /= 0
      self%feces_p = index(self%feces_elements, 'p') /= 0
      
      call self%get_parameter(self%carcass_elements, 'carcass_elements', '', 'elemental composition of carcass e.g. "cnp"')
      self%carcass_c = index(self%carcass_elements, 'c') /= 0
      self%carcass_n = index(self%carcass_elements, 'n') /= 0
      self%carcass_p = index(self%carcass_elements, 'p') /= 0
      
      call self%get_parameter(self%respiration_elements, 'respiration_elements', '', 'elemental composition of respiration e.g. "cnp"')
      self%resp_c = index(self%respiration_elements, 'c') /= 0
      self%resp_n = index(self%respiration_elements, 'n') /= 0
      self%resp_p = index(self%respiration_elements, 'p') /= 0
            
      ! Calculate gww_targetc: gww to target C unit conversion factor 
      self%gww_mmolc = self%gww_gc * 12.01_rk * 1._rk/1000._rk ! Redfield unit conversion from mmol C m-3 to gWW m-3
      
      ! Calculate gww_targetn and gww_targetp: gww to target N and P unit conversion factors
      self%gww_FISHn = self%gww_mmolc / self%qnc                   ! gww to target N unit
      self%gww_FISHp = self%gww_mmolc / self%qpc                   ! gww to target P unit
      
      ! Register model parameters and variables here.
      
      ! convert to per second
      h     = h/seconds_per_year !NIAZ change to scale factor at parameter registration
      gamma = gamma/seconds_per_year !NIAZ change to scale factor at parameter registration
      !kk    = kk/seconds_per_year !NIAZ change to scale factor at parameter registration - also, not used?
      
      !select case (setupidx)
      !   case ()
      !call setupbasic(szprod, lzprod, bprodin, dfbot, depth, Tp, Tb)
            !call setupbasic2(szprod, lzprod, bprodin, dfbot,nStages, depth, Tp, Tb,etaMature,Fmax,etaF,bET_val)
      !   case()
       !call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,3,Tp,Tm,Tb,depth,self%photic_depth,etamature,shelfdepth,visual,Fmax,etaF)
      
      ! Resource forcing is unused in coupled mode; fixed values keep setupVertical2 on a defined branch
      szprod = 100._rk; lzprod = 100._rk; bprodin = 5._rk; dfbot = -1._rk; dfpho = -1._rk
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
      allocate(self%zoo_availability(self%n_zoo))
      self%zoo_availability = 1._rk
      self%n_l = 0
      self%n_s = 0
      do i = 1, self%n_zoo ! find the number of large and small zooplankton models, register dependencies to both types,
          write (i_str,'(i0)') i
          if (self%n_zoo > 1) then !if more than one zoopankton
              call self%get_parameter(self%zoo_size(i), 'zoo_size'//trim(i_str), '', 'small or large zooplankton? "s" for small, "l" for large') ! register whether current model is large or small
              call self%get_parameter(self%zoo_availability(i), 'zoo_'//trim(i_str)//'_availability', '', 'availability of zooplankton prey '//trim(i_str)//' (between 0 and 1)', default=1._rk, minimum=0._rk, maximum=1._rk)
              if (self%zoo_size(i) == 's') then ! if small
                  self%n_s = self%n_s + 1 ! iterate counter
                  !write (ns_str,'(i0)') self%n_s 
                  !call self%register_model_dependency('small_zooplankton_target'//TRIM(ns_str)) ! register model dependency
              else if (self%zoo_size(i) == 'l') then ! if large
                  self%n_l = self%n_l + 1  ! iterate counter
                  !write (nl_str,'(i0)') self%n_l 
                  !call self%register_model_dependency('large_zooplankton_target'//TRIM(nl_str)) ! register model dependency
              else
                  call self%fatal_error('initialize', 'zoo_size'//trim(i_str)//' must be "s" or "l"')
              end if
          else if (self%n_zoo == 1) then ! if only one zooplankton
              !call self%register_model_dependency('zooplankton_target')
              self%n_l = 1
              self%n_s = 1
          end if   
      end do
      
      ! availability per small and per large zooplankton prey, in the order in which they appear in zoo_size<i>
      ! (small_zooplankton_target<j> is the j-th "s" entry, large_zooplankton_target<j> the j-th "l" entry)
      if (self%n_zoo > 1) then
          allocate(self%smzoo_availability(self%n_s), self%lgzoo_availability(self%n_l))
          self%smzoo_availability = pack(self%zoo_availability, self%zoo_size == 's')
          self%lgzoo_availability = pack(self%zoo_availability, self%zoo_size == 'l')
      end if

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
              
              call self%register_mapped_model_dependency(self%id_smzoo_int(i), 'small_zooplankton_target'//trim(i_str), proportional_change=.true., domain=domain_bottom,id_w=self%dummy_w)
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
      
      ! register model dependencies for recipient nutrient and POM/detritus pools
      if (self%resp_c) call self%register_model_dependency("respiration_c_target")
      if (self%resp_n) call self%register_model_dependency("respiration_n_target")
      if (self%resp_p) call self%register_model_dependency("respiration_p_target")
      if (self%feces_c) call self%register_model_dependency("feces_c_target")
      if (self%feces_n) call self%register_model_dependency("feces_n_target")
      if (self%feces_p) call self%register_model_dependency("feces_p_target")
      if (self%carcass_c) call self%register_model_dependency("carcasses_c_target")
      if (self%carcass_n) call self%register_model_dependency("carcasses_n_target")
      if (self%carcass_p) call self%register_model_dependency("carcasses_p_target")

      ! Register the pelagic pools that receive the fish fluxes, one per flux type
      self%flux_on = [self%resp_c, self%resp_n, self%resp_p, self%feces_c, self%feces_n, self%feces_p, self%carcass_c, self%carcass_n, self%carcass_p]
      do k = 1, n_flux
         if (.not. self%flux_on(k)) cycle
         call self%register_state_dependency(self%id_flux_target(k), trim(flux_names(k))//'_pool', 'mmol m-3', trim(flux_names(k))//' receiving pool')
         select case (mod(k-1, 3))
            case (0)
               call self%request_coupling_to_model(self%id_flux_target(k), trim(flux_names(k))//'_target', standard_variables%total_carbon)
            case (1)
               call self%request_coupling_to_model(self%id_flux_target(k), trim(flux_names(k))//'_target', standard_variables%total_nitrogen)
            case (2)
               call self%request_coupling_to_model(self%id_flux_target(k), trim(flux_names(k))//'_target', standard_variables%total_phosphorus)
         end select
         self%id_flux_target(k)%sms%link%target%source = source_do_column  ! sources are added in do_column (see type_projected_rate_distributor in depth_mapping.F90)
      end do
      call self%register_dependency(self%id_z, standard_variables%depth)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)

      ! Determine whether benthos will be modelled internally or externally
      call self%get_parameter(self%internal_benthos, 'internal_benthos', '', 'benthos modelled internally (True) or from external zoobenthos model (False)')
      if (self%internal_benthos) then ! if no coupling to external zoobenthos prey, setup internal zoobenthos variables, and couple to detrital fluxes instead.
          !Register fake deposition pools, for POM/detritus to dump their sedimentation in, and register a dependency on its ROC, to extract the flux
          call self%get_parameter(self%benk, 'benk', 'g m-2', 'benthos carrying capacity')  
          call self%get_parameter(self%ben_eff, 'ben_eff', '-', 'detrital fraction going to benthos ')
          ! deposition into det_pool_* is a direct variable coupling, so FABM applies no unit conversion: these factors convert the depositing model's units to mmol (e.g. 1/12.011 for ERSEM carbon in mg C)
          call self%get_parameter(self%det_c_scale, 'det_c_scale', 'mmol C (deposited unit)-1', 'factor converting deposited carbon to mmol C', default=1._rk)
          call self%get_parameter(self%det_n_scale, 'det_n_scale', 'mmol N (deposited unit)-1', 'factor converting deposited nitrogen to mmol N', default=1._rk)
          call self%get_parameter(self%det_p_scale, 'det_p_scale', 'mmol P (deposited unit)-1', 'factor converting deposited phosphorus to mmol P', default=1._rk)
          call self%get_parameter(self%det_s_scale, 'det_s_scale', 'mmol Si (deposited unit)-1', 'factor converting deposited silicon to mmol Si', default=1._rk)
          
          call self%register_dependency(self%id_det_flux_c,  'det_flux_c',   'mmol m-2 s-1',  'total bottom c deposition flux')
          call self%register_dependency(self%id_det_flux_n,  'det_flux_n',   'mmol m-2 s-1',  'total bottom n deposition flux')
          call self%register_dependency(self%id_det_flux_p,  'det_flux_p',   'mmol m-2 s-1',  'total bottom p deposition flux')
          call self%register_dependency(self%id_det_flux_s,  'det_flux_s',   'mmol m-2 s-1',  'total bottom s deposition flux')
                    
          call self%register_diagnostic_variable(self%id_det_pool_c, 'det_pool_c','mmol m-2', 'target bottom c pool for pelagic tracer deposition',source=source_constant, output=output_none, act_as_state_variable=.true.)
          call self%register_diagnostic_variable(self%id_det_pool_n, 'det_pool_n','mmol m-2', 'target bottom n pool for pelagic tracer deposition',source=source_constant, output=output_none, act_as_state_variable=.true.)
          call self%register_diagnostic_variable(self%id_det_pool_p, 'det_pool_p','mmol m-2', 'target bottom p pool for pelagic tracer deposition',source=source_constant, output=output_none, act_as_state_variable=.true.)
          call self%register_diagnostic_variable(self%id_det_pool_s, 'det_pool_s','mmol m-2', 'target bottom s pool for pelagic tracer deposition',source=source_constant, output=output_none, act_as_state_variable=.true.)
                    
          call self%request_coupling(self%id_det_flux_c,'./det_pool_c_sms_tot')
          call self%request_coupling(self%id_det_flux_n,'./det_pool_n_sms_tot')
          call self%request_coupling(self%id_det_flux_p,'./det_pool_p_sms_tot')
          call self%request_coupling(self%id_det_flux_s,'./det_pool_s_sms_tot')
          
          ! register dependencies on the real bottom detritus pools
          call self%get_parameter(self%bottom_detritus_elements, 'bottom_detritus_elements', '', 'elemental composition of bottom detritus e.g. "cnps"')
          self%bot_det_c = index(self%bottom_detritus_elements, 'c') /= 0
          self%bot_det_n = index(self%bottom_detritus_elements, 'n') /= 0
          self%bot_det_p = index(self%bottom_detritus_elements, 'p') /= 0
          self%bot_det_s = index(self%bottom_detritus_elements, 's') /= 0
          if (self%bot_det_c) then
              call self%register_model_dependency(self%id_det_target_c, 'det_bot_c_target_real') ! register model dependency
              call self%register_state_dependency(self%id_det_bot_c, 'det_bot_c', 'mmol m-2', 'real benthic c pool')
              call self%request_coupling_to_model(self%id_det_bot_c, self%id_det_target_c, standard_variables%total_carbon)
          endif
          if (self%bot_det_n) then
              call self%register_model_dependency(self%id_det_target_n, 'det_bot_n_target_real') ! register model dependency
              call self%register_state_dependency(self%id_det_bot_n, 'det_bot_n', 'mmol m-2', 'real benthic n pool')
              call self%request_coupling_to_model(self%id_det_bot_n, self%id_det_target_n, standard_variables%total_nitrogen)
          endif
          if (self%bot_det_p) then
              call self%register_model_dependency(self%id_det_target_p, 'det_bot_p_target_real') ! register model dependency
              call self%register_state_dependency(self%id_det_bot_p, 'det_bot_p', 'mmol m-2', 'real benthic p pool')
              call self%request_coupling_to_model(self%id_det_bot_p, self%id_det_target_p, standard_variables%total_phosphorus)
          endif
          if (self%bot_det_s) then
              call self%register_model_dependency(self%id_det_target_s, 'det_bot_s_target_real') ! register model dependency
              call self%register_state_dependency(self%id_det_bot_s, 'det_bot_s', 'mmol m-2', 'real benthic si pool')
              call self%request_coupling_to_model(self%id_det_bot_s, self%id_det_target_s, standard_variables%total_silicate)
          endif
          
          ! register dependencies on pelagic nutrient pools recieving benthos excretion flux
          if (self%resp_c) then
              call self%register_state_dependency(self%id_nut_c, 'nut_c', 'mmol m-3', 'pelagic inorganic c pool')
              call self%request_coupling_to_model(self%id_nut_c, "respiration_c_target", standard_variables%total_carbon)
          endif
          if (self%resp_n) then
              call self%register_state_dependency(self%id_nut_n, 'nut_n', 'mmol m-3', 'pelagic inorganic n pool')
              call self%request_coupling_to_model(self%id_nut_n, "respiration_n_target", standard_variables%total_nitrogen)
          endif
          if (self%resp_p) then
              call self%register_state_dependency(self%id_nut_p, 'nut_p', 'mmol m-3', 'pelagic inorganic p pool')
              call self%request_coupling_to_model(self%id_nut_p, "respiration_p_target", standard_variables%total_phosphorus)
          endif
                         
          ! register statevariable for internal zoobenthos 
          call self%register_state_variable(self%id_benthos, 'benthos', 'g m-2', 'benthos biomass', initial_value=bent_ini, minimum=0.0_rk)
          
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
              call self%get_parameter(self%ben_availability(i), 'ben_'//TRIM(i_str)//'_availability', '', 'availability of zoobenthos prey'//TRIM(i_str)//' (between 0 and 1)',minimum=0._rk, maximum=1._rk)
              call self%register_model_dependency(self%id_benprey(i), 'zoobenthos_target'//TRIM(i_str)) ! register model dependency
        
              ! register dependencies on zoobenthos C,N,P
              call self%register_dependency(self%id_ben_c(i), 'zoobenthos_c_'//trim(i_str), 'mmol m-2', 'zoobenthos '//trim(i_str)//' c biomass')
              call self%register_dependency(self%id_ben_n(i), 'zoobenthos_n_'//trim(i_str), 'mmol m-2', 'zoobenthos '//trim(i_str)//' n biomass')
              call self%register_dependency(self%id_ben_p(i), 'zoobenthos_p_'//trim(i_str), 'mmol m-2', 'zoobenthos '//trim(i_str)//' p biomass')
              
              call self%request_coupling_to_model(self%id_ben_c(i), self%id_benprey(i), standard_variables%total_carbon)
              call self%request_coupling_to_model(self%id_ben_n(i), self%id_benprey(i), standard_variables%total_nitrogen)
              call self%request_coupling_to_model(self%id_ben_p(i), self%id_benprey(i), standard_variables%total_phosphorus)
          enddo
      endif

      allocate(self%id_fish(nGrid-nResources))
      allocate(self%id_flux(nGrid-nResources, n_flux), self%id_flux_in(nGrid-nResources, n_flux))

      allocate(self%id_fish_total_biomass(nGroups))
      allocate(self%id_fish_catch(nGroups))
      allocate(self%id_fish_F(nGroups))
      allocate(self%Fsource(nGroups))
      allocate(self%id_Fext(nGroups))
      allocate(self%fft_short_names(nGrid-nResources))
      allocate(self%size_numstr(nGrid-nResources))

      ! Prepare pre-computed matrices for depth-dependent fish vertical distribution and predation preference (MUST be done before creating child models)
      call initialize_theta_matrices(self, szprod, lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, etamature, visual, Fmax, etaF)  ! Read from file or compute theta/vertical distribution matrices

      !assign fish background mortality to FEISTY
      !mort0(idxF:nGrid) = self%bgmort/seconds_per_year !NIAZ change to scale factor at parameter registration
      !mortF(idxF:nGrid) = mortF/seconds_per_year !NIAZ change to scale factor at parameter registration

      ! registration of fish state variables and diagnostic variables
      n_names = 1
      do i = 1, nGroups
         write (i_str,'(i0)') i
         ! read functional type longnames from yaml
         !call self%get_parameter(fft_long_name,  'fft_long_name'//trim(i_str),  units='', long_name='fish functional type long name',  default="fish_functional_type_"//trim(i_str)) !NIAZ: not used
         ! read functional type shortnames from yaml
         call self%get_parameter(fft_short_name, 'fft_short_name_'//trim(i_str), units='', long_name='fish functional type short name', default="fft_"//trim(i_str))
         ! per functional type fishing; defaults to the global Fmax/etaF. Must come after initialize_theta_matrices,
         ! since setupVertical2 (called there) resets fishing mortality of all groups to the global values.
         call self%get_parameter(self%Fsource(i), 'Fsource_'//trim(i_str), '', 'source of fishing mortality for '//trim(fft_short_name)// &
              ' (0: Fmax parameter, 1: external maximum fishing mortality F_'//trim(fft_short_name)//' with FEISTY selectivity)', default=0)
         call self%get_parameter(etaF_i, 'etaF_'//trim(i_str), '-', 'fraction of maximum size with 50% fishing selectivity for '//trim(fft_short_name), default=etaF)
         select case (self%Fsource(i))
         case (0)
            call self%get_parameter(Fmax_i, 'Fmax_'//trim(i_str), 'yr-1', 'maximum fishing mortality of '//trim(fft_short_name), default=Fmax)
            call setFishing(Fmax_i, etaF_i, [i])
         case (1)
            ! store the selectivity curve only (Fmax = 1); do_bottom scales it by F_<fft> every time step
            call setFishing(1._rk, etaF_i, [i])
            call self%register_dependency(self%id_Fext(i), 'F_'//trim(fft_short_name), 'yr-1', 'maximum fishing mortality of '//trim(fft_short_name))
         case default
            call self%fatal_error('initialize', 'Fsource_'//trim(i_str)//' must be 0 or 1')
         end select
         ! read the initial value of each functional type from yaml
         !call self%get_parameter(self%b_ini, 'b_ini_'//trim(i_str), 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  

         do j = ixStart(i)-nResources, ixEnd(i)-nResources
               write (j_str,'(i0)') j
               write (size_number_str,'(i0)') j+nResources-ixStart(i)+1
               ! assign the initial biomass to each fish size class of a functional type.
               call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(size_number_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(size_number_str)//'_biomass', initial_value=self%b_ini, minimum=0.0_rk)
               
               self%fft_short_names(n_names) = fft_short_name
               self%size_numstr(n_names) = size_number_str
               n_names = n_names + 1
         end do
            
         call self%register_diagnostic_variable(self%id_fish_total_biomass(i), trim(fft_short_name)//'_totB', 'g m-2', trim(fft_short_name)//'_total_biomass')
         call self%register_diagnostic_variable(self%id_fish_catch(i), trim(fft_short_name)//'_catch', 'g m-2 yr-1', trim(fft_short_name)//'_catch (wet weight)')
         call self%register_diagnostic_variable(self%id_fish_F(i), trim(fft_short_name)//'_F', 'yr-1', 'fishing mortality of the largest size class of '//trim(fft_short_name))
      end do 
      
      !sets up biogeochemical coupling and vertical distribution for each fish size class
      do i = 1, nGrid-nResources !NIAZ change this to a nested loop over functional types and size classes like in the previous loop, and use the functional type short name for the variable and coupling names
         write (i_str,'(i0)') i
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_fish(i), scale_factor = 1._rk/self%gww_mmolc)     !if (self%nutrient_in_c) NIAZ: remove the target scaling and use redfield.
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fish(i), scale_factor = 1._rk/self%gww_FISHn)   !if (self%nutrient_in_n) NIAZ: remove the target scaling and use redfield.
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_fish(i), scale_factor = 1._rk/self%gww_FISHp) !if (self%nutrient_in_p) NIAZ: remove the target scaling and use redfield.
         
         ! Fluxes to the BGC model: stored as bottom diagnostics in do_bottom, distributed over the water column in do_column
         do k = 1, n_flux
            if (.not. self%flux_on(k)) cycle
            flux_name = trim(flux_names(k))//'_fish_'//trim(self%fft_short_names(i))//'_size'//trim(self%size_numstr(i))
            call self%register_diagnostic_variable(self%id_flux(i,k), trim(flux_name), 'mmol m-2 s-1', trim(flux_name), source=source_do_bottom, output=output_none)
            call self%register_dependency(self%id_flux_in(i,k), trim(flux_name)//'_in', 'mmol m-2 s-1', trim(flux_name))
            call self%request_coupling(self%id_flux_in(i,k), './'//trim(flux_name))
         end do
         
      end do
       
       !temperature and depth
       call self%register_dependency(self%id_temp, standard_variables%temperature)
       call self%register_dependency(self%id_temp_vertmean_100m, vertical_mean(self%id_temp, maximum_depth=100._rk))
       call self%register_dependency(self%id_temp_vertmean_500to1500m, vertical_mean(self%id_temp, minimum_depth=500._rk, maximum_depth=1500._rk))
       call self%register_diagnostic_variable(self%id_temp_vertmean_100m_diag, 'temp_vertmean_100m', 'degree_C', 'vertical mean temperature above 100 m')
       call self%register_dependency(self%id_bottom_depth, standard_variables%bottom_depth)
       
!   -----declared in setup.f90-----
      allocate (respiration(nGrid))
      allocate (carcasses(nGrid))
      allocate (feces(nGrid))
      
   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_feisty_fabm), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: p1
      integer  :: iGroup, i, k, istate, iprey, closest_bottom_depth_idx
      real(rk),dimension(nGrid-nResources)       :: fish, pred_morts, w_sz, w_lz, w_ben
      real(rk),dimension(nGrid-nResources,n_flux) :: flux
      real(rk)                                   :: growth_frac
      real(rk)                                   :: F_ext
      real(rk)                                   :: benthos_g_gww
      real(rk),dimension(nGrid)                  :: uin, dudt
      real(rk)                                   :: temp_100m,temp_bottom, temp_500to1500m, bottom_depth
      real(rk)                                   :: photic_depth_local  ! local photic depth adjusted per grid cell
      real(rk),dimension(self%n_s)               :: smzoo_c, smzoo_n, smzoo_p, smzoo_ww_rf
      real(rk),dimension(self%n_l)               :: lgzoo_c, lgzoo_n, lgzoo_p, lgzoo_ww_rf
      real(rk),dimension(1)                      :: zoo_c, zoo_n, zoo_p 
      real(rk),dimension(self%n_ben)             :: benthos_c, benthos_n, benthos_p, ben_ww_rf
      real(rk)                                   :: smzoo_ww_tot, lgzoo_ww_tot, smzoo_mort_ww, lgzoo_mort_ww
      real(rk)                                   :: excess_c_sz, excess_n_sz, excess_p_sz, excess_c_lz, excess_n_lz, excess_p_lz, ben_ww_tot, excess_c_ben, excess_n_ben, excess_p_ben, ben_mort_ww
      real(rk)                                   :: det_flux_c, det_flux_n, det_flux_p, det_flux_s, det_flux_ww, excess_c_det, excess_n_det, excess_p_det, det_surplus_c, det_surplus_n, det_surplus_p
      
      _BOTTOM_LOOP_BEGIN_

         ! Get depth-integrated predator biomass       
         do i = 1, nGrid-nResources
         _GET_BOTTOM_(self%id_fish(i), fish(i))
         end do      
         !print*,(sum(fish))
         
         ! get depth-integrated zooplankton
         if (self%n_zoo > 1) then ! if theres more than one zooplankton, get all biomasses and then sum up for each size group.
             do i = 1, self%n_s
                _GET_BOTTOM_(self%id_smzoo_c(i), smzoo_c(i)) ! NIAZ: These are in mmol m-2
                _GET_BOTTOM_(self%id_smzoo_n(i), smzoo_n(i)) ! NIAZ: These are in mmol m-2
                _GET_BOTTOM_(self%id_smzoo_p(i), smzoo_p(i)) ! NIAZ: These are in mmol m-2
             enddo
                          
             ! calculate WW of limiting nutrient for each small zooplankton, scaled by its availability
             smzoo_ww_rf = merge(min(merge(max(smzoo_c,0._rk)*self%smzoo_availability*self%gww_mmolc, huge(0._rk), smzoo_c /= 0._rk), &
                                     merge(max(smzoo_n,0._rk)*self%smzoo_availability*self%gww_FISHn, huge(0._rk), smzoo_n /= 0._rk), &
                                     merge(max(smzoo_p,0._rk)*self%smzoo_availability*self%gww_FISHp, huge(0._rk), smzoo_p /= 0._rk)), 0._rk, smzoo_c /= 0._rk .or. smzoo_n /= 0._rk .or. smzoo_p /= 0._rk)
             do i = 1, self%n_l
                _GET_BOTTOM_(self%id_lgzoo_c(i), lgzoo_c(i)) ! NIAZ: These are in mmol m-2
                _GET_BOTTOM_(self%id_lgzoo_n(i), lgzoo_n(i)) ! NIAZ: These are in mmol m-2
                _GET_BOTTOM_(self%id_lgzoo_p(i), lgzoo_p(i)) ! NIAZ: These are in mmol m-2
             enddo
                          
             ! calculate WW of limiting nutrient for each large zooplankton, scaled by its availability
             lgzoo_ww_rf = merge(min(merge(max(lgzoo_c,0._rk)*self%lgzoo_availability*self%gww_mmolc, huge(0._rk), lgzoo_c /= 0._rk), &
                                     merge(max(lgzoo_n,0._rk)*self%lgzoo_availability*self%gww_FISHn, huge(0._rk), lgzoo_n /= 0._rk), &
                                     merge(max(lgzoo_p,0._rk)*self%lgzoo_availability*self%gww_FISHp, huge(0._rk), lgzoo_p /= 0._rk)), 0._rk, lgzoo_c /= 0._rk .or. lgzoo_n /= 0._rk .or. lgzoo_p /= 0._rk)
             
             
         else ! if theres only one registered zooplankton prey, get its biomass, and then split evenly between small and large. NIAZ: Ask Ken about this!!!!
             _GET_BOTTOM_(self%id_zoo_c, zoo_c(1)) ! NIAZ: These are in mmol m-2
             _GET_BOTTOM_(self%id_zoo_n, zoo_n(1)) ! NIAZ: These are in mmol m-2
             _GET_BOTTOM_(self%id_zoo_p, zoo_p(1)) ! NIAZ: These are in mmol m-2
             
             !calculate WW of limiting nutrient
             smzoo_ww_rf = merge(min(merge(max(zoo_c,0._rk)*self%gww_mmolc, huge(0._rk), zoo_c /= 0._rk), &
                                     merge(max(zoo_n,0._rk)*self%gww_FISHn, huge(0._rk), zoo_n /= 0._rk), &
                                     merge(max(zoo_p,0._rk)*self%gww_FISHp, huge(0._rk), zoo_p /= 0._rk)), 0._rk, zoo_c /= 0._rk .or. zoo_n /= 0._rk .or. zoo_p /= 0._rk) / 2._rk
             lgzoo_ww_rf = smzoo_ww_rf
             
             
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
             ben_ww_rf = merge(min(merge(max(benthos_c,0._rk)*self%ben_availability*self%gww_mmolc, huge(0._rk), benthos_c /= 0._rk), &
                                   merge(max(benthos_n,0._rk)*self%ben_availability*self%gww_FISHn, huge(0._rk), benthos_n /= 0._rk), &
                                   merge(max(benthos_p,0._rk)*self%ben_availability*self%gww_FISHp, huge(0._rk), benthos_p /= 0._rk)), 0._rk, benthos_c /= 0._rk .or. benthos_n /= 0._rk .or. benthos_p /= 0._rk)
             ben_ww_tot = sum(ben_ww_rf)
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
         if (self%depth_array(closest_bottom_depth_idx) .lt. (photic_depth_local + 500._dp)) dvm = self%depth_array(closest_bottom_depth_idx) 
         if (self%depth_array(closest_bottom_depth_idx) .le. self%shelfdepth) dvm = 0._dp

         !update temperature effect
         call updateTempV2(temp_100m, temp_500to1500m, temp_bottom, dvm, bottom_depth, photic_depth_local, ixmedium, ixlarge)
         do iGroup = 1, nGroups
             group(iGroup)%spec%V=group(iGroup)%spec%Vsave*fTempV(ixStart(iGroup):ixEnd(iGroup))
             group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmaxsave*fTempV(ixStart(iGroup):ixEnd(iGroup))
             group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolismsave*fTempmV(ixStart(iGroup):ixEnd(iGroup))
         end do
         call set2vec
         call scale_mortality_rates(self)

         ! external fishing mortality (Fsource=1): group%spec%mortF holds the selectivity curve, scale it by F_<fft>.
         ! Same operation order as the parameter path, (psiF*F)/seconds_per_year, so a constant input reproduces Fsource=0 exactly.
         do iGroup = 1, nGroups
            if (self%Fsource(iGroup) == 1) then
               _GET_BOTTOM_(self%id_Fext(iGroup), F_ext)
               mortF(ixStart(iGroup):ixEnd(iGroup)) = (group(iGroup)%spec%mortF*max(F_ext,0._rk))/seconds_per_year
            end if
         end do

         ! FEISTY derivatives
         uin= [smzoo_ww_tot,lgzoo_ww_tot ,ben_ww_tot,0._rk ,fish] ! NIAZ: at this point, zooplankton and benthos biomasses (in WW) must be calculated and aggregated.
         call calcderivatives(uin, dudt)
         
         smzoo_ww_tot = max(smzoo_ww_tot,tiny(0.0_rk))
         lgzoo_ww_tot = max(lgzoo_ww_tot,tiny(0.0_rk))
         ben_ww_tot = max(ben_ww_tot,tiny(0.0_rk))
         !print*,dudt
         

         !do i = 1, nGrid
         uin = max(0._rk , uin)
         !end do
         !! recalculate each predator’s share of the grazing on each prey, to distribute the excess nutrient excretion. NIAZ: nGrid-nResources
         pred_morts = Cmax(nResources+1:nGrid) * V(nResources+1:nGrid) / (Enc(nResources+1:nGrid) + Cmax(nResources+1:nGrid)) * uin(nResources+1:nGrid)
         
         w_sz  = theta(nResources+1:nGrid,1) * pred_morts
         w_lz  = theta(nResources+1:nGrid,2) * pred_morts
         w_ben = theta(nResources+1:nGrid,3) * pred_morts
         w_sz  = w_sz  / max(sum(w_sz),  tiny(0._rk))
         w_lz  = w_lz  / max(sum(w_lz),  tiny(0._rk))
         w_ben = w_ben / max(sum(w_ben), tiny(0._rk))
                  
         ! NIAZ: Zooplankton removal should be scaled to the relative WW of the most and least limiting prey nutrients - the difference is then excreted.
         smzoo_mort_ww = dudt(1)/smzoo_ww_tot
         lgzoo_mort_ww = dudt(2)/lgzoo_ww_tot
         
         !calculate total zooplankton mortality as dudt(1)/smzoo_ww_tot, then apply mortality to state variables of prey
         if (self%n_zoo > 1) then ! if multiple zooplankton prey types
             ! update small zooplankton state var ROCs
             do iprey = 1, self%n_s
                do istate = 1, size(self%id_smzoo_int(iprey)%bottom_state)
                   _GET_BOTTOM_(self%id_smzoo_int(iprey)%bottom_state(istate), p1)
                   _ADD_BOTTOM_SOURCE_(self%id_smzoo_int(iprey)%bottom_state(istate), smzoo_mort_ww * self%smzoo_availability(iprey) * p1) ! NIAZ: this assumes dudt(1) is already negative.
                end do
             end do
             ! update large zooplankton state var ROCs
             do iprey = 1, self%n_l
                do istate = 1, size(self%id_lgzoo_int(iprey)%bottom_state)
                   _GET_BOTTOM_(self%id_lgzoo_int(iprey)%bottom_state(istate), p1)
                   _ADD_BOTTOM_SOURCE_(self%id_lgzoo_int(iprey)%bottom_state(istate), lgzoo_mort_ww * self%lgzoo_availability(iprey) * p1) ! NIAZ: this assumes dudt(2) is already negative.
                end do
             end do
             
             ! calculate excretion/respiration due to stochiometry
             excess_c_sz = sum(smzoo_c*self%smzoo_availability*self%gww_mmolc * (-smzoo_mort_ww) + dudt(1)*smzoo_ww_rf/smzoo_ww_tot) / self%gww_mmolc ! gww biomass * specific loss rate - loss
             excess_n_sz = sum(smzoo_n*self%smzoo_availability*self%gww_FISHn * (-smzoo_mort_ww) + dudt(1)*smzoo_ww_rf/smzoo_ww_tot) / self%gww_FISHn ! gww biomass * specific loss rate - loss
             excess_p_sz = sum(smzoo_p*self%smzoo_availability*self%gww_FISHp * (-smzoo_mort_ww) + dudt(1)*smzoo_ww_rf/smzoo_ww_tot) / self%gww_FISHp ! gww biomass * specific loss rate - loss

             excess_c_lz = sum(lgzoo_c*self%lgzoo_availability*self%gww_mmolc * (-lgzoo_mort_ww) + dudt(2)*lgzoo_ww_rf/lgzoo_ww_tot) / self%gww_mmolc ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-2
             excess_n_lz = sum(lgzoo_n*self%lgzoo_availability*self%gww_FISHn * (-lgzoo_mort_ww) + dudt(2)*lgzoo_ww_rf/lgzoo_ww_tot) / self%gww_FISHn ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-2
             excess_p_lz = sum(lgzoo_p*self%lgzoo_availability*self%gww_FISHp * (-lgzoo_mort_ww) + dudt(2)*lgzoo_ww_rf/lgzoo_ww_tot) / self%gww_FISHp ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-2
         else ! if only one zooplankton prey type
             ! Apply mortality to prey state variables
             do istate = 1, size(self%id_zoo_int%bottom_state)
                _GET_BOTTOM_(self%id_zoo_int%bottom_state(istate), p1)
                _ADD_BOTTOM_SOURCE_(self%id_zoo_int%bottom_state(istate), (lgzoo_mort_ww + smzoo_mort_ww)/2._rk * p1) ! NIAZ: this assumes dudt(X) is already negative.
             end do
             
             ! calculate excretion/respiration due to stochiometry
             excess_c_sz = sum(zoo_c*self%gww_mmolc * (-smzoo_mort_ww)/2._rk + dudt(1)) / self%gww_mmolc
             excess_n_sz = sum(zoo_n*self%gww_FISHn * (-smzoo_mort_ww)/2._rk + dudt(1)) / self%gww_FISHn
             excess_p_sz = sum(zoo_p*self%gww_FISHp * (-smzoo_mort_ww)/2._rk + dudt(1)) / self%gww_FISHp
             excess_c_lz = sum(zoo_c*self%gww_mmolc * (-lgzoo_mort_ww)/2._rk + dudt(2)) / self%gww_mmolc
             excess_n_lz = sum(zoo_n*self%gww_FISHn * (-lgzoo_mort_ww)/2._rk + dudt(2)) / self%gww_FISHn
             excess_p_lz = sum(zoo_p*self%gww_FISHp * (-lgzoo_mort_ww)/2._rk + dudt(2)) / self%gww_FISHp
         endif
      
         !calculate total zoobenthos mortality, then apply mortality to state variables of prey
         ben_mort_ww = dudt(3)/ben_ww_tot
         if (self%internal_benthos) then ! if coupling to detritus fluxes, running an internal benthos model.
             ! get the deposition fluxes
             _GET_BOTTOM_(self%id_det_flux_c, det_flux_c)
             _GET_BOTTOM_(self%id_det_flux_n, det_flux_n)
             _GET_BOTTOM_(self%id_det_flux_p, det_flux_p)
             _GET_BOTTOM_(self%id_det_flux_s, det_flux_s)
             ! convert to mmol m-2 s-1; everything below (incl. the surplus returned via the standard-variable det_bot_* couplings) is then in mmol
             det_flux_c = det_flux_c*self%det_c_scale
             det_flux_n = det_flux_n*self%det_n_scale
             det_flux_p = det_flux_p*self%det_p_scale
             det_flux_s = det_flux_s*self%det_s_scale

             ! calculate WW of limiting nutrient for detritus flux           
             det_flux_ww = merge(min(merge(max(det_flux_c,0._rk)*self%gww_mmolc, huge(0._rk), det_flux_c /= 0._rk), &
                                     merge(max(det_flux_n,0._rk)*self%gww_FISHn, huge(0._rk), det_flux_n /= 0._rk), &
                                     merge(max(det_flux_p,0._rk)*self%gww_FISHp, huge(0._rk), det_flux_p /= 0._rk)), 0._rk, det_flux_c /= 0._rk .or. det_flux_n /= 0._rk .or. det_flux_p /= 0._rk)
             
             ! calculate zoobenthos growth 
             growth_frac = self%ben_eff*max(0._rk, 1._rk - ben_ww_tot/self%benk)
             benthos_g_gww = growth_frac*det_flux_ww
             
             ! update zoobenthos ROC
             _ADD_BOTTOM_SOURCE_(self%id_benthos, benthos_g_gww + dudt(3))
             
             ! calcualate excess nutrients, in case detritus stochiometry differs from RF ratios
             excess_c_det = growth_frac*det_flux_c - benthos_g_gww/self%gww_mmolc
             excess_n_det = growth_frac*det_flux_n - benthos_g_gww/self%gww_FISHn
             excess_p_det = growth_frac*det_flux_p - benthos_g_gww/self%gww_FISHp
             
             det_surplus_c = (1._rk - growth_frac)*det_flux_c
             det_surplus_n = (1._rk - growth_frac)*det_flux_n
             det_surplus_p = (1._rk - growth_frac)*det_flux_p
             
             ! update bottom detritus
             if (self%bot_det_c) then
                 _ADD_BOTTOM_SOURCE_(self%id_det_bot_c, det_surplus_c)
             endif
             if (self%bot_det_n) then
                 _ADD_BOTTOM_SOURCE_(self%id_det_bot_n, det_surplus_n)
             endif
             if (self%bot_det_p) then
                 _ADD_BOTTOM_SOURCE_(self%id_det_bot_p, det_surplus_p)
             endif
             if (self%bot_det_s) then
                 _ADD_BOTTOM_SOURCE_(self%id_det_bot_s, det_flux_s)
             endif
             
             ! update pelagic nutrient pools
             if (self%resp_c) then
                 _ADD_BOTTOM_FLUX_(self%id_nut_c, excess_c_det)
             endif
             if (self%resp_n) then
                 _ADD_BOTTOM_FLUX_(self%id_nut_n, excess_n_det)
             endif
             if (self%resp_p) then
                 _ADD_BOTTOM_FLUX_(self%id_nut_p, excess_p_det)
             endif
             
             ! no excess nutrients from fish predation, since internal benthos keep RF ratios
             excess_c_ben = 0.0_rk
             excess_n_ben = 0.0_rk 
             excess_p_ben = 0.0_rk
         else ! if using external benthos model(s)
             do iprey = 1, self%n_ben
                 do istate = 1, size(self%id_benprey(iprey)%bottom_state)
                     _GET_BOTTOM_(self%id_benprey(iprey)%bottom_state(istate), p1)
                     _ADD_BOTTOM_SOURCE_(self%id_benprey(iprey)%bottom_state(istate), ben_mort_ww * self%ben_availability(iprey) * p1) ! NIAZ: this assumes dudt(3) is already negative.
                 end do
             end do
                 
             ! calculate excretion/respiration due to stochiometry
             excess_c_ben = sum(benthos_c*self%ben_availability*self%gww_mmolc * (-ben_mort_ww) + dudt(3)*ben_ww_rf/ben_ww_tot) / self%gww_mmolc ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-2
             excess_n_ben = sum(benthos_n*self%ben_availability*self%gww_FISHn * (-ben_mort_ww) + dudt(3)*ben_ww_rf/ben_ww_tot) / self%gww_FISHn ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-2
             excess_p_ben = sum(benthos_p*self%ben_availability*self%gww_FISHp * (-ben_mort_ww) + dudt(3)*ben_ww_rf/ben_ww_tot) / self%gww_FISHp ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-2
         endif
         
         do i = 1, nGrid-nResources
             ! fish dynamics
             _ADD_BOTTOM_SOURCE_(self%id_fish(i), dudt(i+nResources) )!gww/m2
             
             ! fluxes to the BGC model (mmol m-2 s-1); do_column distributes them over the water column
             flux(i, flux_resp_c)    = respiration(i+nResources) / self%gww_mmolc + excess_c_sz*w_sz(i) + excess_c_lz*w_lz(i) + excess_c_ben*w_ben(i)
             flux(i, flux_resp_n)    = respiration(i+nResources) / self%gww_FISHn + excess_n_sz*w_sz(i) + excess_n_lz*w_lz(i) + excess_n_ben*w_ben(i)
             flux(i, flux_resp_p)    = respiration(i+nResources) / self%gww_FISHp + excess_p_sz*w_sz(i) + excess_p_lz*w_lz(i) + excess_p_ben*w_ben(i)
             flux(i, flux_feces_c)   = feces(i+nResources) / self%gww_mmolc      !gww/m2 to target C/m2
             flux(i, flux_feces_n)   = feces(i+nResources) / self%gww_FISHn      !gww/m2 to mmol N/m2
             flux(i, flux_feces_p)   = feces(i+nResources) / self%gww_FISHp      !gww/m2 to mmol P/m2
             flux(i, flux_carcass_c) = carcasses(i+nResources) / self%gww_mmolc  !gww/m2 to target C/m2
             flux(i, flux_carcass_n) = carcasses(i+nResources) / self%gww_FISHn  !gww/m2 to mmol N/m2
             flux(i, flux_carcass_p) = carcasses(i+nResources) / self%gww_FISHp  !gww/m2 to mmol P/m2
             do k = 1, n_flux
                if (self%flux_on(k)) then
                   _SET_BOTTOM_DIAGNOSTIC_(self%id_flux(i,k), flux(i,k))
                end if
             end do
         end do
         
         do i=1,nGroups
            _SET_BOTTOM_DIAGNOSTIC_(self%id_fish_total_biomass(i), totBiomass(i))
            ! catch = fishing mortality x biomass (as in R calcYield); uin is clamped >= 0 above, mortF is in s-1
            _SET_BOTTOM_DIAGNOSTIC_(self%id_fish_catch(i), sum(mortF(ixStart(i):ixEnd(i))*uin(ixStart(i):ixEnd(i)))*seconds_per_year)
            _SET_BOTTOM_DIAGNOSTIC_(self%id_fish_F(i), mortF(ixEnd(i))*seconds_per_year)
         end do
      _BOTTOM_LOOP_END_
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      ! Distributes the fluxes of all fish size classes over the water column, following the FEISTY vertical distributions
      class (type_feisty_fabm), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: bottom_depth, z_center, h, scale_factor, max_depth
      real(rk),dimension(nGrid-nResources)        :: w, w_int
      real(rk),dimension(nGrid-nResources,n_flux) :: flux
      integer  :: i, k, closest_bottom_depth_idx

      ! Get the fluxes of each fish size class (mmol m-2 s-1), calculated in do_bottom
      flux = 0.0_rk
      do k = 1, n_flux
         if (.not. self%flux_on(k)) cycle
         do i = 1, nGrid-nResources
            _GET_HORIZONTAL_(self%id_flux_in(i,k), flux(i,k))
         end do
      end do

      ! Find the closest pre-computed depth profile, and the scale factor that compresses it into the water column
      _GET_HORIZONTAL_(self%id_bottom_depth, bottom_depth)
      closest_bottom_depth_idx = minloc(abs(self%depth_array - bottom_depth), dim=1)
      scale_factor = bottom_depth / self%depth_array(closest_bottom_depth_idx)
      max_depth = min(self%depth_array(closest_bottom_depth_idx), real(size(self%vertical_distribution_matrix_centers, 1), rk))

      ! First pass: depth-integrated weights of each fish size class (m), so that each flux is conserved exactly
      w_int = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_z, z_center)
         _GET_(self%id_h, h)
         call feisty_layer_weights(self%vertical_distribution_matrix_centers(:,:,closest_bottom_depth_idx), z_center, h, scale_factor, max_depth, w)
         w_int = w_int + w * h
      _DOWNWARD_LOOP_END_
      where (w_int <= 0.0_rk) w_int = huge(w_int)  ! no profile in this column: its flux is dropped (should not happen)

      ! Second pass: add the summed fluxes to the pelagic pools (mmol m-3 s-1)
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_z, z_center)
         _GET_(self%id_h, h)
         call feisty_layer_weights(self%vertical_distribution_matrix_centers(:,:,closest_bottom_depth_idx), z_center, h, scale_factor, max_depth, w)
         do k = 1, n_flux
            if (self%flux_on(k)) then
               _ADD_SOURCE_(self%id_flux_target(k), sum(flux(:,k) * w / w_int))
            end if
         end do
      _DOWNWARD_LOOP_END_
   end subroutine do_column

   subroutine initialize_theta_matrices(self, szprod, lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, etamature, visual, Fmax, etaF)
      class (type_feisty_fabm), intent(inout) :: self
      real(rk), intent(in) :: szprod, lzprod, bprodin, dfbot, dfpho
      integer, intent(in)  :: nStages
      real(rk), intent(in) :: Tp, Tm, Tb, etamature, visual, Fmax, etaF
      logical :: file_exists
      integer :: dims(3)

      call allocate_vertical_arrays(self) 
      self%theta_matrix = 0._rk
      self%vertical_distribution_matrix_centers = 0._rk ! NIAZ: changed from -999._rk on advice from Claude

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
      integer :: i_depth
    
      if (self%n_depth > int(self%max_depth)) call self%fatal_error('initialize', 'n_depth must not exceed max_depth (m)') ! NIAZ: Added to avoid infinite looping below
      
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
      !! Verify all depths are unique (quality check) NIAZ: These could never trigger with the above logic
      !is_unique = .true.
      !do i_depth = 1, self%n_depth - 1
      !   do j_depth = i_depth + 1, self%n_depth
      !      if (self%depth_array(i_depth) == self%depth_array(j_depth)) then
      !         print*, "WARNING: depth_array still contains duplicates after removal loop!"
      !         is_unique = .false.
      !         exit
      !      end if
      !   end do
      !   if (.not. is_unique) exit
      !end do
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
         
         call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,nStages,Tp,Tm,Tb,depth_local,photic_depth_local,etamature,self%shelfdepth,visual,Fmax,etaF)
         
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
      !if (self%bgmort .ne. mort0(nGrid)) then
         mort0(idxF:nGrid) = self%bgmort!/seconds_per_year 
      !else
         !mort0(idxF:nGrid) = mort0(idxF:nGrid)/seconds_per_year 
      !end if
      mortF(idxF:nGrid) = mortF(idxF:nGrid)/seconds_per_year 
   end subroutine scale_mortality_rates

   ! Add model subroutines here.

end module
