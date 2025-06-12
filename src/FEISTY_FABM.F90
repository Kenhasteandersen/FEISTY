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
      type (type_bottom_state_variable_id),         allocatable :: id_smpel(:), id_mesopel(:), id_lgpel(:), id_midp(:), id_dem(:)
      ! The vertical distributios of the fish:
      type (type_vertical_distribution_id),         allocatable :: id_smpel_w(:), id_mesopel_w(:), id_lgpel_w(:), id_midp_w(:), id_dem_w(:)
      ! The benthos state variable:
      type (type_bottom_state_variable_id)                      :: id_benthos

      ! Dependency IDs for the state variables in the biogeochemical model:
      ! Small zooplankton carbon (c), nitrogen (n), and phosphorus (p):
      ! (note that we register as a bottom variable, but it is actually summed over the water column)
      type (type_bottom_dependency_id)                         :: id_smzoo_c, id_smzoo_n, id_smzoo_p
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
 
      ! Coupling; pointer to which model contains the zooplankton state varaiable
      ! that we will integrater over the water column:
      type (type_model_id)                                      :: id_zooplankton 
      
   contains
      procedure :: initialize
      !  procedure :: do_bottom
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
      integer :: i
      character(len=10)  :: strindex
      
      class (type_vertical_depth_range), pointer :: depth_distribution
      
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
      
      call setupbasic(szprod, lzprod, bprodin, dfbot, depth, Tp, Tb)
      
      ! Register state variables
      !allocate(self%id_u(nGrid))
      !call self%register_state_variable(self%id_u(1), 'u1', 'g m-2', 'biomass', initial_value=smz_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(2), 'u2', 'g m-2', 'biomass', initial_value=lgz_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(3), 'u3', 'g m-2', 'biomass', initial_value=smbent_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(4), 'u4', 'g m-2', 'biomass', initial_value=lgbent_ini, minimum=0.0_rk)      

      allocate(self%id_smpel(2))
      allocate(self%id_smpel_w(2))! allpcate size of vertical distribution of small pelagics
      do i = 1, 2 !small pelagic
         write (strindex,'(i0)') i
         call self%register_state_variable(self%id_smpel(i), 'smpel'//trim(strindex), 'g m-2', 'small pelagic biomass', initial_value=b_ini, minimum=0.0_rk)
         call self%set_variable_property(self%id_smpel(i), 'disable_transport', .true.)
         
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_smpel(i), scale_factor = gwwC)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_smpel(i), scale_factor = gwwC*NC)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_smpel(i), scale_factor = gwwC*PC)
         
         call self%register_vertical_distribution(self%id_smpel_w(i),'smpel_w'//trim(strindex))
         allocate(depth_distribution)
         call self%add_child(depth_distribution, 'habitat_smpel_w'//trim(strindex))
         call self%request_coupling(self%id_smpel_w(i), 'habitat_smpel_w'//trim(strindex)//'/w')
      end do

      allocate(self%id_lgpel(3))
      allocate(self%id_lgpel_w(3))
      do i = 1, 3 !large pelagic
         write (strindex,'(i0)') i
         call self%register_state_variable(self%id_lgpel(i), 'lgpel'//trim(strindex), 'g m-2', 'large pelagic biomass', initial_value=b_ini, minimum=0.0_rk)
         call self%set_variable_property(self%id_lgpel(i), 'disable_transport', .true.)
         
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_lgpel(i), scale_factor = gwwC)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_lgpel(i), scale_factor = gwwC*NC)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_lgpel(i), scale_factor = gwwC*PC)

         call self%register_vertical_distribution(self%id_lgpel_w(i),'lgpel_w'//trim(strindex))
         allocate(depth_distribution)
         call self%add_child(depth_distribution, 'habitat_lgpel_w'//trim(strindex))
         call self%request_coupling(self%id_lgpel_w(i), 'habitat_lgpel_w'//trim(strindex)//'/w')
      
      end do
      
      allocate(self%id_dem(3))
      allocate(self%id_dem_w(3))
      do i = 1, 3 !demersal fish
         write (strindex,'(i0)') i
         call self%register_state_variable(self%id_dem(i), 'dem'//trim(strindex), 'g m-2', 'demersal biomass', initial_value=b_ini, minimum=0.0_rk)
         call self%set_variable_property(self%id_dem(i), 'disable_transport', .true.)
         
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_dem(i), scale_factor = gwwC)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_dem(i), scale_factor = gwwC*NC)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_dem(i), scale_factor = gwwC*PC)

         call self%register_vertical_distribution(self%id_dem_w(i),'dem_w'//trim(strindex))
         allocate(depth_distribution)         
         call self%add_child(depth_distribution, 'habitat_dem_w'//trim(strindex))
         call self%request_coupling(self%id_dem_w(i), 'habitat_dem_w'//trim(strindex)//'/w')      
      end do
      
      
      !benthos register
      call self%register_state_variable(self%id_benthos,'benthos', 'g m-2', 'benthos biomass', initial_value=smbent_ini, minimum=0.0_rk)
      
       ! Depth-averaged dependencies
       !call self%register_dependency(self%id_temp, 'temp', 'degrees_Celsius', 'depthaveraged temperature')
       call self%register_dependency(self%id_smzoo_c, 'smzoo_c', 'mmol C m-2', 'depth-integrated small zooplankton carbon')
       call self%register_dependency(self%id_smzoo_n, 'smzoo_n', 'mmol N m-2', 'depth-integrated small zooplankton nitrogen')
       call self%register_dependency(self%id_smzoo_p, 'smzoo_p', 'mmol P m-2', 'depth-integrated small zooplankton phosphorus') 
       
       
       call self%register_state_dependency(self%id_excre_n,'excretion_N', 'mmol N m-2', 'excretion nitrogen')
       call self%register_state_dependency(self%id_excre_p,'excretion_p', 'mmol C m-2', 'excretion phosphorus')
       
       call self%register_state_dependency(self%id_respiration_c,'respiration_c', 'mmol C m-2', 'respiration carbon')
       
       call self%register_state_dependency(self%id_feces_c,'feces_c', 'mmol C m-2', 'feces carbon')
       call self%register_state_dependency(self%id_feces_n,'feces_n', 'mmol N m-2', 'feces nitrogen')
       call self%register_state_dependency(self%id_feces_p,'feces_p', 'mmol P m-2', 'feces phosphorus')

       call self%register_state_dependency(self%id_carcasses_c,'carcasses_c', 'mmol C m-2', 'carcasses carbon')
       call self%register_state_dependency(self%id_carcasses_n,'carcasses_n', 'mmol N m-2', 'carcasses nitrogen')
       call self%register_state_dependency(self%id_carcasses_p,'carcasses_p', 'mmol P m-2', 'carcasses phosphorus')
      
       
       call self%request_mapped_coupling_to_model(self%id_smzoo_c, 'small_zooplankton',standard_variables%total_carbon)
       call self%request_mapped_coupling_to_model(self%id_smzoo_n, 'small_zooplankton',standard_variables%total_nitrogen)
       call self%request_mapped_coupling_to_model(self%id_smzoo_p, 'small_zooplankton',standard_variables%total_phosphorus)

       call self%request_mapped_coupling_to_model(self%id_excre_n, 'excretion',standard_variables%total_nitrogen)
       call self%request_mapped_coupling_to_model(self%id_excre_p, 'excretion',standard_variables%total_phosphorus)  
       
       call self%request_mapped_coupling_to_model(self%id_feces_c, 'feces',standard_variables%total_carbon)
       call self%request_mapped_coupling_to_model(self%id_feces_n, 'feces',standard_variables%total_nitrogen)
       call self%request_mapped_coupling_to_model(self%id_feces_p, 'feces',standard_variables%total_phosphorus)
       
       call self%request_mapped_coupling_to_model(self%id_carcasses_c, 'carcasses',standard_variables%total_carbon)
       call self%request_mapped_coupling_to_model(self%id_carcasses_n, 'carcasses',standard_variables%total_nitrogen)
       call self%request_mapped_coupling_to_model(self%id_carcasses_p, 'carcasses',standard_variables%total_phosphorus)
       
       call self%register_mapped_model_dependency(self%id_zooplankton, 'small_zooplankton', proportional_change=.true., domain=domain_bottom)
       
       
       
      
!      do i = idxF, nGrid
!         write (strindex,'(i0)') i
!         call self%register_state_variable(self%id_u(i), 'u'//trim(strindex), 'g m-2', 'biomass', initial_value=b_ini, minimum=0.0_rk)
!               call self%set_variable_property(self%id_u(i), 'disable_transport', .true.)
!      end do
      
   end subroutine initialize

   ! Add model subroutines here.

end module
