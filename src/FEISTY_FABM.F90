#include "fabm_driver.h"

module FEISTY_FABM

   use globals
   use spectrum
   use fish
   use setup
   use fabm_types
   use fabm_particle
   use fabm_builtin_depth_mapping   
   
   implicit none

   private

   type, extends(type_depth_integrated_particle), public :: type_feisty_fabm
      ! State variables for fish (small pelagics, mesopelagics, large pelagics, 
      ! mid-water predators, and demersal fish):
      type (type_bottom_state_variable_id),         allocatable :: id_fish(:)!, id_smpel(:), id_mesopel(:), id_lgpel(:), id_midp(:), id_dem(:)
      ! The vertical distributios of the fish:
      type (type_vertical_distribution_id),         allocatable :: id_fish_w(:)!, id_smpel_w(:), id_mesopel_w(:), id_lgpel_w(:), id_midp_w(:), id_dem_w(:)
      ! The benthos state variable:
      type (type_bottom_state_variable_id)                      :: id_benthos

      ! Dependency IDs for the state variables in the biogeochemical model:
      ! Small zooplankton carbon (c), nitrogen (n), and phosphorus (p):
      ! (note that we register as a bottom variable, but it is actually summed over the water column)
      type (type_bottom_dependency_id),         allocatable   :: id_smzoo_fish_c(:), id_smzoo_fish_n(:), id_smzoo_fish_p(:)
      !type (type_bottom_dependency_id),         allocatable   :: id_smzoo_smpel_c(:), id_smzoo_smpel_n(:), id_smzoo_smpel_p(:)
      !type (type_bottom_dependency_id),         allocatable   :: id_smzoo_lgpel_c(:), id_smzoo_lgpel_n(:), id_smzoo_lgpel_p(:)      
      !type (type_bottom_dependency_id),         allocatable   :: id_smzoo_dem_c(:), id_smzoo_dem_n(:), id_smzoo_dem_p(:)      
      ! type (type_bottom_dependency_id)                         :: id_lgzoo_c,id_lgzoo_n,id_lgzoo_p
      ! type (type_bottom_dependency_id)                           :: id_temp

      ! Dependency IDs for the outputs from FEISTY
      ! Excretion of nitrogen (n) and phosphorus (p):
      type (type_bottom_state_variable_id)                      :: id_excre_n, id_excre_p
      ! Respiration (CO2 / DIC):
      type (type_bottom_state_variable_id)                      :: id_respiration_c
      ! Fecal pellets (for detritus / POM):
      type (type_bottom_state_variable_id)                      :: id_feces_c, id_feces_n, id_feces_p
      ! Carcasses:
      type (type_bottom_state_variable_id)                      :: id_carcasses_c, id_carcasses_n, id_carcasses_p 
      
      !type (type_bottom_dependency_id),    allocatable      :: id_excre_smpel_n(:), id_excre_smpel_p(:), id_respiration_smpel_c(:),id_feces_smpel_c(:), id_feces_smpel_n(:), id_feces_smpel_p(:),id_carcasses_smpel_c(:), id_carcasses_smpel_n(:), id_carcasses_smpel_p(:)
      !type (type_bottom_dependency_id),    allocatable      :: id_excre_lgpel_n(:), id_excre_lgpel_p(:), id_respiration_lgpel_c(:),id_feces_lgpel_c(:), id_feces_lgpel_n(:), id_feces_lgpel_p(:),id_carcasses_lgpel_c(:), id_carcasses_lgpel_n(:), id_carcasses_lgpel_p(:) 
      !type (type_bottom_dependency_id),    allocatable      :: id_excre_dem_n(:), id_excre_dem_p(:), id_respiration_dem_c(:),id_feces_dem_c(:), id_feces_dem_n(:), id_feces_dem_p(:),id_carcasses_dem_c(:), id_carcasses_dem_n(:), id_carcasses_dem_p(:) 
      type (type_bottom_dependency_id),    allocatable      :: id_excre_fish_n(:), id_excre_fish_p(:), id_respiration_fish_c(:),id_feces_fish_c(:), id_feces_fish_n(:), id_feces_fish_p(:),id_carcasses_fish_c(:), id_carcasses_fish_n(:), id_carcasses_fish_p(:) 
      ! Coupling; pointer to which model contains the zooplankton state varaiable
      ! that we will integrater over the water column:
      !type (type_model_id)                                      :: id_zooplankton 
      
   contains
      procedure :: initialize
      procedure :: do_bottom
      ! Reference model procedures here.
   end type type_feisty_fabm
   
    ! Redfieldian N:C and P:C ratios of predator biomass
    real(rk), parameter :: NC = 16.0_rk/ 106.0_rk
    real(rk), parameter :: PC = 1.0_rk / 106.0_rk
    !from gww to C
    real(rk), parameter :: gwwC = 1._rk/9._rk

contains

   subroutine initialize(self, configunit)
      class (type_feisty_fabm), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
      
      real(rk)           :: smz_ini, lgz_ini, smbent_ini, lgbent_ini, b_ini
      real(rk)           :: szprod, lzprod, bprodin, dfbot, depth, Tp, Tb
      integer :: i,j
      character(len=100)  :: strindex, i_str, j_str, fft_long_name,fft_short_name
      
      class (type_vertical_depth_range), pointer :: depth_distribution
      
      !call self%type_depth_integrated_particle%initialize(configunit)
      
      ! Register model parameters  
      ! Fish pysilogical parameters. See input.nml in the original FEISTY R package and Fortran library
      call self%get_parameter(h, 'h', 'g^nn yr-1', 'Max. consumption coefficient', default=20._rk)  
      call self%get_parameter(nn, 'nn', '-', 'Metabolic exponent', default=-0.25_rk)  
      call self%get_parameter(gamma, 'gamma', 'm2 g^q yr-1', 'Coef. for clearance rate', default=70._rk)  
      call self%get_parameter(q, 'q', '-', 'Clearance rate exponent', default=-0.2_rk)  
      call self%get_parameter(kk, 'kk', 'g^p yr-1', 'Metabolism coefficient', default=4._rk)  
      call self%get_parameter(p, 'p', '-', 'Metabolism exponent', default=-0.175_rk)  
      call self%get_parameter(epsAssim, 'epsAssim', '-', 'Assimilation efficiency', default=0.7_rk)  
      call self%get_parameter(epsRepro, 'epsRepro', '-', 'Reproduction & recruitment efficiency', default=0.01_rk)  
      
      call self%get_parameter(beta, 'beta', '-', 'Beta parameter for size-based predation preference', default=400._rk)  
      call self%get_parameter(sigma, 'sigma', '-', 'Sigma parameter for size-based predation preference', default=1.3_rk)  
      call self%get_parameter(mMin, 'mMin', 'g', 'Minimum fish mass (boundary of the grid)', default=0.001_rk)  
      call self%get_parameter(mMedium, 'mMedium', 'g', 'Medium fish central mass for feeding preference calc', default=0.5_rk)  
      call self%get_parameter(mLarge, 'mLarge', 'g', 'Large fish central mass for feeding preference calc', default=250._rk) 
      
      call self%get_parameter(lbenk, 'lbenk', 'g m-2', 'Large benthos carry capacity', default=0._rk)  
      call self%get_parameter(szoog, 'szoog', 'yr-1', 'Small zooplankton growth rate', default=1._rk)  
      call self%get_parameter(lzoog, 'lzoog', 'yr-1', 'Large zooplankton growth rate', default=1._rk)  
      call self%get_parameter(sbeng, 'sbeng', 'yr-1', 'Small benthos growth rate', default=1._rk)  
      call self%get_parameter(lbeng, 'lbeng', 'yr-1', 'Large benthos growth rate', default=0._rk)  
      ! predation preference coefficient
      call self%get_parameter(thetaS, 'thetaS', '-', 'Medium fish preference for small zooplankton', default=0.25_rk)  
      call self%get_parameter(thetaA, 'thetaA', '-', 'Large fish preference for medium forage fish', default=0.5_rk)  
      call self%get_parameter(thetaD, 'thetaD', '-', 'Preference of large demersal on pelagic prey', default=0.75_rk)  
      ! get initial values
      call self%get_parameter(smz_ini,'smz_ini', 'g m-2',     'initial small mesozooplankton biomass', default=1.e2_rk)  
      call self%get_parameter(lgz_ini, 'lgz_ini', 'g m-2', 'initial large mesozooplankton biomass', default=1.e2_rk)  
      call self%get_parameter(smbent_ini, 'smbent_ini', 'g m-2', 'initial small benthos biomass', default=5._rk)  
      call self%get_parameter(lgbent_ini, 'lgbent_ini', 'g m-2', 'initial large benthos biomass', default=0._rk)        
      call self%get_parameter(b_ini, 'b_ini', 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  
      ! get input values
      call self%get_parameter(szprod, 'szprod', 'g m-2 yr-1', 'small mesozooplankton production', default=100._rk)  
      call self%get_parameter(lzprod, 'lzprod', 'g m-2 yr-1', 'large mesozooplankton production', default=100._rk)  
      call self%get_parameter(bprodin, 'bprodin', 'g m-2 yr-1', 'benthos production', default=5._rk)  
      call self%get_parameter(dfbot, 'dfbot', 'g m-2 yr-1', 'detrital flux reaching the bottom', default=-1._rk)  
      call self%get_parameter(depth, 'depth', 'm', 'water column depth', default=100._rk)  
      call self%get_parameter(Tp, 'Tp', 'Celsius', 'pelagic layer averaged temperature', default=10._rk)  
      call self%get_parameter(Tb, 'Tb', 'Celsius', 'bottom layer depth temperature', default=8._rk)      
      ! Register model parameters and variables here.
      
      !select case (setupidx)
      !   case ()
      call setupbasic(szprod, lzprod, bprodin, dfbot, depth, Tp, Tb)
      !   case()
      !      
      !end select      
      
      
      ! Register state variables
      !allocate(self%id_u(nGrid))
      !call self%register_state_variable(self%id_u(1), 'u1', 'g m-2', 'biomass', initial_value=smz_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(2), 'u2', 'g m-2', 'biomass', initial_value=lgz_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(3), 'u3', 'g m-2', 'biomass', initial_value=smbent_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(4), 'u4', 'g m-2', 'biomass', initial_value=lgbent_ini, minimum=0.0_rk)      

      call self%register_model_dependency("small_zooplankton")
      call self%register_model_dependency("excretion")
      call self%register_model_dependency("respiration")
      call self%register_model_dependency("feces")
      call self%register_model_dependency("carcasses")
      
      allocate(self%id_fish(nGrid-nResources))
      allocate(self%id_fish_w(nGrid-nResources))! allpcate size of vertical distribution of small pelagics     
      allocate(self%id_smzoo_fish_c(nGrid-nResources), self%id_smzoo_fish_n(nGrid-nResources), self%id_smzoo_fish_p(nGrid-nResources))
      allocate(self%id_excre_fish_n(nGrid-nResources), self%id_excre_fish_p(nGrid-nResources))
      allocate(self%id_respiration_fish_c(nGrid-nResources))
      allocate(self%id_feces_fish_c(nGrid-nResources), self%id_feces_fish_n(nGrid-nResources), self%id_feces_fish_p(nGrid-nResources))
      allocate(self%id_carcasses_fish_c(nGrid-nResources), self%id_carcasses_fish_n(nGrid-nResources), self%id_carcasses_fish_p(nGrid-nResources))
      
      do i = 1, nGroups
         write (i_str,'(i0)') i
         ! read functional type longnames from yaml
         call self%get_parameter(fft_long_name,  'fft_long_name'//trim(i_str),  units='', long_name='fish functional type long name',  default="fish_functional_type_"//trim(i_str))
         ! read functional type shortnames from yaml
         call self%get_parameter(fft_short_name, 'fft_short_name_'//trim(i_str), units='', long_name='fish functional type short name', default="fft_"//trim(i_str))
         ! read the initial value of each functional type from yaml
         call self%get_parameter(b_ini, 'b_ini_'//trim(i_str), 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  

            do j = ixStart(i)-nResources, ixEnd(i)-nResources
               write (j_str,'(i0)') j
               ! assign the initial biomass to each fish size class of a functional type.
                call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(j_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(i_str)//'_biomass', initial_value=b_ini, minimum=0.0_rk)

            end do
      end do 
      
      do i = 1, nGrid-nResources
         write (i_str,'(i0)') i
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_fish(i), scale_factor = gwwC)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fish(i), scale_factor = gwwC*NC)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_fish(i), scale_factor = gwwC*PC)
         
         call self%register_vertical_distribution(self%id_fish_w(i),'fish_w'//trim(i_str))
         allocate(depth_distribution)
         call self%add_child(depth_distribution, 'habitat_fish_w'//trim(i_str))
         call self%request_coupling(self%id_fish_w(i), 'habitat_fish_w'//trim(i_str)//'/'//'w')
         
         call self%register_dependency(self%id_smzoo_fish_c(i), 'smzoo_c', 'mmol C m-2', 'depth-integrated small zooplankton carbon')
         call self%register_dependency(self%id_smzoo_fish_n(i), 'smzoo_n', 'mmol N m-2', 'depth-integrated small zooplankton nitrogen')
         call self%register_dependency(self%id_smzoo_fish_p(i), 'smzoo_p', 'mmol P m-2', 'depth-integrated small zooplankton phosphorus')   
       
         call self%register_dependency(self%id_excre_fish_n(i),'excretion_n', 'mmol N m-2', 'excretion nitrogen')
         call self%register_dependency(self%id_excre_fish_p(i),'excretion_p', 'mmol C m-2', 'excretion phosphorus')
       
         call self%register_dependency(self%id_respiration_fish_c(i),'respiration_c', 'mmol C m-2', 'respiration carbon')
       
         call self%register_dependency(self%id_feces_fish_c(i),'feces_c', 'mmol C m-2', 'feces carbon')
         call self%register_dependency(self%id_feces_fish_n(i),'feces_n', 'mmol N m-2', 'feces nitrogen')
         call self%register_dependency(self%id_feces_fish_p(i),'feces_p', 'mmol P m-2', 'feces phosphorus')
       
         call self%register_dependency(self%id_carcasses_fish_c(i),'carcasses_c', 'mmol C m-2', 'carcasses carbon')
         call self%register_dependency(self%id_carcasses_fish_n(i),'carcasses_n', 'mmol N m-2', 'carcasses nitrogen')
         call self%register_dependency(self%id_carcasses_fish_p(i),'carcasses_p', 'mmol P m-2', 'carcasses phosphorus')
         
         call self%request_mapped_coupling_to_model(self%id_smzoo_fish_c(i), 'small_zooplankton_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_smzoo_fish_n(i), 'small_zooplankton_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_smzoo_fish_p(i), 'small_zooplankton_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
         call self%couplings%set_string('small_zooplankton_fish_'//trim(i_str), "small_zooplankton")
       
         call self%request_mapped_coupling_to_model(self%id_excre_fish_n(i), 'excretion_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_excre_fish_p(i), 'excretion_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))  
         call self%couplings%set('excretion_fish_'//trim(i_str), "excretion")

         call self%request_mapped_coupling_to_model(self%id_respiration_fish_c(i),'respiration_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         call self%couplings%set_string('respiration_fish_'//trim(i_str), "respiration")

         call self%request_mapped_coupling_to_model(self%id_feces_fish_c(i), 'feces_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_feces_fish_n(i), 'feces_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_feces_fish_p(i), 'feces_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
         call self%couplings%set_string('feces_fish_'//trim(i_str), "feces")

         call self%request_mapped_coupling_to_model(self%id_carcasses_fish_c(i), 'carcasses_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_carcasses_fish_n(i), 'carcasses_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_carcasses_fish_p(i), 'carcasses_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
         call self%couplings%set_string('carcasses_fish_'//trim(i_str), "carcasses")         
         
      end do
      
      
      
      !benthos register
      call self%register_state_variable(self%id_benthos, 'benthos', 'g m-2', 'benthos biomass', initial_value=smbent_ini, minimum=0.0_rk)
      
       ! Depth-averaged dependencies
       !call self%register_dependency(self%id_temp, 'temp', 'degrees_Celsius', 'depthaveraged temperature')
       !call self%register_dependency(self%id_smzoo_c, 'smzoo_c', 'mmol C m-2', 'depth-integrated small zooplankton carbon')
       !call self%register_dependency(self%id_smzoo_n, 'smzoo_n', 'mmol N m-2', 'depth-integrated small zooplankton nitrogen')
       !call self%register_dependency(self%id_smzoo_p, 'smzoo_p', 'mmol P m-2', 'depth-integrated small zooplankton phosphorus') 
       
       
       !call self%register_state_dependency(self%id_excre_n,'excretion_N', 'mmol N m-2', 'excretion nitrogen')
       !call self%register_state_dependency(self%id_excre_p,'excretion_p', 'mmol C m-2', 'excretion phosphorus')
       !
       !call self%register_state_dependency(self%id_respiration_c,'respiration_c', 'mmol C m-2', 'respiration carbon')
       !
       !call self%register_state_dependency(self%id_feces_c,'feces_c', 'mmol C m-2', 'feces carbon')
       !call self%register_state_dependency(self%id_feces_n,'feces_n', 'mmol N m-2', 'feces nitrogen')
       !call self%register_state_dependency(self%id_feces_p,'feces_p', 'mmol P m-2', 'feces phosphorus')
       !
       !call self%register_state_dependency(self%id_carcasses_c,'carcasses_c', 'mmol C m-2', 'carcasses carbon')
       !call self%register_state_dependency(self%id_carcasses_n,'carcasses_n', 'mmol N m-2', 'carcasses nitrogen')
       !call self%register_state_dependency(self%id_carcasses_p,'carcasses_p', 'mmol P m-2', 'carcasses phosphorus')
      
       
       !call self%request_mapped_coupling_to_model(self%id_smzoo_c, 'small_zooplankton',standard_variables%total_carbon, vertical_distribution=self%id_lgpel_w(i))
       !call self%request_mapped_coupling_to_model(self%id_smzoo_n, 'small_zooplankton',standard_variables%total_nitrogen, vertical_distribution=self%id_lgpel_w(i))
       !call self%request_mapped_coupling_to_model(self%id_smzoo_p, 'small_zooplankton',standard_variables%total_phosphorus, vertical_distribution=self%id_lgpel_w(i))
       !
       !call self%request_mapped_coupling_to_model(self%id_excre_n, 'excretion',standard_variables%total_nitrogen)
       !call self%request_mapped_coupling_to_model(self%id_excre_p, 'excretion',standard_variables%total_phosphorus)  
       !
       !call self%request_mapped_coupling_to_model(self%id_respiration_c,'respiration',standard_variables%total_carbon)
       !
       !call self%request_mapped_coupling_to_model(self%id_feces_c, 'feces',standard_variables%total_carbon)
       !call self%request_mapped_coupling_to_model(self%id_feces_n, 'feces',standard_variables%total_nitrogen)
       !call self%request_mapped_coupling_to_model(self%id_feces_p, 'feces',standard_variables%total_phosphorus)
       !
       !call self%request_mapped_coupling_to_model(self%id_carcasses_c, 'carcasses',standard_variables%total_carbon)
       !call self%request_mapped_coupling_to_model(self%id_carcasses_n, 'carcasses',standard_variables%total_nitrogen)
       !call self%request_mapped_coupling_to_model(self%id_carcasses_p, 'carcasses',standard_variables%total_phosphorus)
       
       !call self%register_mapped_model_dependency(self%id_zooplankton, 'small_zooplankton', proportional_change=.true., domain=domain_bottom)
       
       
       
      
!      do i = idxF, nGrid
!         write (strindex,'(i0)') i
!         call self%register_state_variable(self%id_u(i), 'u'//trim(strindex), 'g m-2', 'biomass', initial_value=b_ini, minimum=0.0_rk)
!               call self%set_variable_property(self%id_u(i), 'disable_transport', .true.)
!      end do
      
   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_feisty_fabm), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c, temp, prey_c, prey_n, prey_p, prey_s, w_int
      real(rk) :: ingestion_c, ingestion_n, ingestion_p, prey_loss_rate, p, net_growth
      integer  :: iGroup, i, istate
      real(rk),dimension(ixEnd(1)-ixStart(1)+1)  :: smpel
      real(rk),dimension(ixEnd(2)-ixStart(2)+1)  :: lgpel
      real(rk),dimension(ixEnd(3)-ixStart(3)+1)  :: dem
      real(rk),dimension(nGrid-nResources)       :: fish
      real(rk)                                   :: zoo1, zoo2, benthos1, benthos2
      real(rk),dimension(nGrid)                  :: uin, dudt

      _BOTTOM_LOOP_BEGIN_
         ! Get depth-integrated predator biomass       
         do i = 1, nGrid-nResources
         _GET_BOTTOM_(self%id_fish(i), fish(i))
         end do      
         
         _GET_BOTTOM_(self%id_smzoo_fish_c(1),zoo1)
         _GET_BOTTOM_(self%id_smzoo_fish_c(1),zoo2)
         _GET_BOTTOM_(self%id_smzoo_fish_c(1),benthos1)
         _GET_BOTTOM_(self%id_smzoo_fish_c(1),benthos2)
         !print*,(fish(1))
         print*,(zoo1)
         ! Depth-averaged environmental dependencies and prey concentrations
         !_GET_BOTTOM_(self%id_temp, temp)
         !_GET_BOTTOM_(self%id_prey_c, prey_c)
         !_GET_BOTTOM_(self%id_prey_n, prey_n)
         !_GET_BOTTOM_(self%id_prey_p, prey_p)
!         _GET_BOTTOM_(self%id_smzoo_c, c)
         
         uin= [zoo1,zoo2,benthos1,benthos2,fish]
         call calcderivatives(uin, dudt)
         !_ADD_BOTTOM_SOURCE_(self%id_fish, dudt(idxF:nGrid))

         ! Calculate ingested fluxes of different chemical elements
         ! Predator population growth will be based on the most limiting of these
         !ingestion_c = self%clearance_rate * c * prey_c
         !ingestion_n = self%clearance_rate * c * prey_n
         !ingestion_p = self%clearance_rate * c * prey_p
         !net_growth = min(ingestion_c, ingestion_n / NC, ingestion_p / PC) - self%mortality * c

         ! The specific loss rate of prey is the depth-integrated ingestion,
         ! divided by depth-integrated prey biomass, e.g., ingestion_c / prey_c_int.
         ! In turn, prey_c_int is related to depth-averaged prey as prey_c = prey_c_int / w_int,
         ! with w_int representing the depth-integral weights of the predator's vertical distibution.
         ! Thus, the specific loss rate is ingestion_c / (prey_c * w_int), which simplifies to
         ! clearance_rate * c / w_int (see expression for ingestion_c above)
         !_GET_BOTTOM_(self%id_w%integral, w_int)
         !prey_loss_rate = self%clearance_rate * c / w_int

        ! small pelagic
         !do i = 1,(ixEnd(1)-ixStart(1)+1)
         !_GET_BOTTOM_(self%id_smpel_w(i)%integral, w_int)
         !end do
         
         ! Source term for predator
         !_ADD_BOTTOM_SOURCE_(self%id_c, net_growth)

         ! Apply the same specific loss rate of all state variables of the prey
         !do istate = 1, size(self%id_zooplankton%bottom_state)
         !   _GET_BOTTOM_(self%id_zooplankton%bottom_state(istate), p)
            !_ADD_BOTTOM_SOURCE_(self%id_prey_int%bottom_state(istate), -prey_loss_rate * p)
         !end do

         ! Send unused ingested matter and dead biomass to waste pools
         !_ADD_BOTTOM_SOURCE_(self%id_waste_c, ingestion_c - net_growth)
         !_ADD_BOTTOM_SOURCE_(self%id_waste_n, ingestion_n - net_growth * NC)
         !_ADD_BOTTOM_SOURCE_(self%id_waste_p, ingestion_p - net_growth * PC)

         ! Save diagnostics
         !_SET_BOTTOM_DIAGNOSTIC_(self%id_net_growth, net_growth * 86400.0_rk)
         !_SET_BOTTOM_DIAGNOSTIC_(self%id_prey_loss_rate, prey_loss_rate * 86400.0_rk)
         
      _BOTTOM_LOOP_END_
   end subroutine   
   
   ! Add model subroutines here.

end module
