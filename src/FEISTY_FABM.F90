#include "fabm_driver.h"

module FEISTY_FABM

   use globals                            ! FEISTY
   use spectrum                           ! FEISTY
   use fish                               ! FEISTY
   use setup                              ! FEISTY
   use fabm_types                         ! FABM
   use fabm_particle                      ! FABM
   use fabm_builtin_depth_mapping         ! FABM
   use feisty_fabm_vertical_distribution  ! FEISTY-controlled vertical distributions
   use fabm_expressions                   ! FABM
   
   implicit none

   private

   type, extends(type_depth_integrated_particle), public :: type_feisty_fabm
      ! State variables for fish (small pelagics, mesopelagics, large pelagics, 
      ! mid-water predators, and demersal fish):
      type (type_bottom_state_variable_id),         allocatable :: id_fish(:)
      ! The vertical distributios of the fish:
      type (type_vertical_distribution_id),         allocatable :: id_fish_w(:)!, id_smpel_w(:), id_mesopel_w(:), id_lgpel_w(:), id_midp_w(:), id_dem_w(:)
      type (type_vertical_distribution_id)                      :: dummy_w
      ! The benthos state variable:
      type (type_bottom_state_variable_id)                      :: id_benthos
      type (type_state_variable_id)                             :: id_det,id_nut

      ! Dependency IDs for the state variables in the biogeochemical model:
      ! Small zooplankton carbon (c), nitrogen (n), and phosphorus (p):
      ! (note that we register as a bottom variable, but it is actually summed over the water column)
      !type (type_bottom_state_variable_id),         allocatable   :: id_smzoo_fish_c(:), id_smzoo_fish_n(:), id_smzoo_fish_p(:)
      !type (type_bottom_state_variable_id),         allocatable   :: id_lgzoo_fish_c(:), id_lgzoo_fish_n(:), id_lgzoo_fish_p(:) 

      ! Dependency IDs for the outputs from FEISTY
      type (type_bottom_dependency_id)          :: id_smzoo_c, id_smzoo_n, id_smzoo_p
      type (type_bottom_dependency_id)          :: id_lgzoo_c, id_lgzoo_n, id_lgzoo_p
      ! Excretion of nitrogen (n) and phosphorus (p):
      type (type_bottom_state_variable_id),    allocatable      :: id_excre_fish_n(:), id_excre_fish_p(:)
      ! Respiration (CO2 / DIC):
      type (type_bottom_state_variable_id),    allocatable      :: id_respiration_fish_c(:),id_respiration_fish_n(:),id_respiration_fish_p(:)
      ! Fecal pellets (for detritus / POM):
      type (type_bottom_state_variable_id),    allocatable      :: id_feces_fish_c(:), id_feces_fish_n(:), id_feces_fish_p(:)
      ! Carcasses:
      type (type_bottom_state_variable_id),    allocatable      :: id_carcasses_fish_c(:), id_carcasses_fish_n(:), id_carcasses_fish_p(:)
      ! Coupling; pointer to which model contains the zooplankton state varaiable
      ! that we will integrater over the water column:
      type (type_model_id)                                      :: id_smzoo_int, id_lgzoo_int
      
      type (type_dependency_id)                                 :: id_temp,id_kpar
      type (type_horizontal_dependency_id)                      :: id_temp_vertmean_100m,id_temp_vertmean_500to1500m,id_bottom_depth
      type (type_bottom_diagnostic_variable_id)                 :: id_temp_vertmean_100m_diag
      type (type_bottom_diagnostic_variable_id), allocatable    :: id_fish_total_biomass(:)
      
      real (rk)                              :: b_ini
      real (rk) ,    allocatable             :: depth_array(:), photic_depth_array(:), theta_matrix(:,:,:,:) !(predator,prey,dvm_depth,bottom_depth,photic_depth)
      real (rk) ,    allocatable             :: vertical_distribution_matrix(:,:,:,:) !(depth_level,nFGrid,n_depth,n_photic_depth)
      
      real(rk)                               :: shelfdepth
      real(rk)                               :: max_depth, depth_interval, max_photic_depth, photic_depth_interval
      integer                                :: n_depth, n_photic_depth
      
   contains
      procedure :: initialize
      procedure :: do_bottom
      ! Reference model procedures here.
   end type type_feisty_fabm
   
    ! Redfieldian N:C and P:C ratios of predator biomass
    real(rk), parameter :: CN = 16.0_rk/ 106.0_rk
    real(rk), parameter :: CP = 1.0_rk / 106.0_rk
    !from gww to C
    real(rk), parameter :: gwwC = 1._rk/9._rk
    real(rk), parameter :: gww_mmolC = 1._rk/9_rk /(12.01_rk/1000._rk) ! gww to gwwC to mol C to mmol C
    real(rk), parameter :: gww_mmolN = 1._rk/9_rk /(12.01_rk/1000._rk) *16._rk/106._rk ! gww to mmol C to mmol N
    real(rk), parameter :: gww_mmolP = 1._rk/9_rk /(12.01_rk/1000._rk) /106._rk ! gww to mmol C to mmol P    

contains

   subroutine initialize(self, configunit)
      class (type_feisty_fabm), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
      
      real(rk)           :: smz_ini, lgz_ini, smbent_ini, lgbent_ini
      real(rk)           :: szprod, lzprod, bprodin, dfbot, depth, Tp, Tb
      real(rk)           :: bgmort
      real(rk)           :: dfpho, Tm, photic, etamature, visual, Fmax, etaF!, ssigma, tau, shelfdepth
      integer            :: nStages, bET_val
      integer            :: i,j
      character(len=100) :: strindex, i_str, j_str, size_number_str, fft_long_name,fft_short_name
      logical            :: file_exists
      integer            :: dims(4)
      
      class (type_vertical_depth_range), pointer :: depth_distribution
      class (type_feisty_vertical_distribution), pointer :: feisty_vertical_distribution
      
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
      call self%get_parameter(bgmort, 'bgmort', 'yr-1', 'Fish background mortality', default=0.1_rk) 
      
      call self%get_parameter(beta, 'beta', '-', 'Beta parameter for size-based predation preference', default=400._rk)  
      call self%get_parameter(sigma, 'sigma', '-', 'Sigma parameter for size-based predation preference', default=1.3_rk)  
      call self%get_parameter(mMin, 'mMin', 'g', 'Minimum fish mass (boundary of the grid)', default=0.001_rk)  
      call self%get_parameter(mMedium, 'mMedium', 'g', 'Medium fish central mass for feeding preference calc', default=0.5_rk)   ! for vertical distribution 10
      call self%get_parameter(mLarge, 'mLarge', 'g', 'Large fish central mass for feeding preference calc', default=250._rk)     ! for vertical distribution 5000
      
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
      call self%get_parameter(photic, 'photic', 'm', 'photic zone depth', default=150._rk)
      call self%get_parameter(etamature, 'etamature', '-', 'the coefficient determines the fish size with a 50% maturity level', default=0.25_rk)
      call self%get_parameter(self%shelfdepth, 'shelfdepth', 'm', 'continental shelf depth', default=250._rk)
      call self%get_parameter(visual, 'visual', '-', 'the coefficient determines the visual ability of fish', default=1.5_rk) ! 1.5:visual predator or 1:non-visual predator. Be careful to use other values.
      call self%get_parameter(Fmax, 'Fmax', 'yr-1', 'Maximum fishing mortality', default=0._rk)
      call self%get_parameter(etaF, 'etaF', '-', 'the coefficient determines the fish size with 50% fishing selectivity', default=0.05_rk)
      call self%get_parameter(ssigma, 'ssigma', '-', 'width of initial vertical distribution', default=10._rk)
      call self%get_parameter(tau, 'tau', '-', 'increase in width', default=10._rk)
      
      call self%get_parameter(bET_val, 'bET_val', '-', 'whether turn on the effective temperature (integer 1 or 0)', default=1)
      
      call self%get_parameter(self%max_depth, 'max_depth', 'm', 'maximum depth of the model domain', default=6000._rk)
      call self%get_parameter(self%depth_interval, 'depth_interval', 'm', 'depth interval for vertical distribution', default=50._rk)
      call self%get_parameter(self%max_photic_depth, 'max_photic_depth', 'm', 'maximum photic depth of the model domain', default=200._rk)
      call self%get_parameter(self%photic_depth_interval, 'photic_depth_interval', 'm', 'photic depth interval for vertical distribution', default=10._rk)! 10 meters?
      self%n_depth        = int(self%max_depth/self%depth_interval)
      self%n_photic_depth = int(self%max_photic_depth/self%photic_depth_interval)
      
      ! Register model parameters and variables here.
      
      ! convert to per second
      h     = h/365._rk/86400._rk
      gamma = gamma/365._rk/86400._rk
      kk    = kk/365._rk/86400._rk
      
      !select case (setupidx)
      !   case ()
      !call setupbasic(szprod, lzprod, bprodin, dfbot, depth, Tp, Tb)
            !call setupbasic2(szprod, lzprod, bprodin, dfbot,nStages, depth, Tp, Tb,etaMature,Fmax,etaF,bET_val)
      !   case()
       !call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,3,Tp,Tm,Tb,depth,photic,etamature,shelfdepth,visual,Fmax,etaF)
       call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,9,Tp,Tm,Tb,depth,photic,etamature,self%shelfdepth,visual,Fmax,etaF)
      !      
      
      !end select     
      
      !assign fish background mortality to FEISTY
      mort0(idxF:nGrid) = bgmort/365._rk/86400._rk
      mortF = mortF/365._rk/86400._rk
      
      ! Register state variables
      !allocate(self%id_u(nGrid))
      !call self%register_state_variable(self%id_u(1), 'u1', 'g m-2', 'biomass', initial_value=smz_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(2), 'u2', 'g m-2', 'biomass', initial_value=lgz_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(3), 'u3', 'g m-2', 'biomass', initial_value=smbent_ini, minimum=0.0_rk)
      !call self%register_state_variable(self%id_u(4), 'u4', 'g m-2', 'biomass', initial_value=lgbent_ini, minimum=0.0_rk)      

      call self%register_model_dependency("small_zooplankton")
      call self%register_model_dependency("large_zooplankton")
      call self%register_model_dependency("excretion")
      call self%register_model_dependency("respiration")
      call self%register_model_dependency("feces")
      call self%register_model_dependency("carcasses")
      
      allocate(self%id_fish(nGrid-nResources))
      allocate(self%id_fish_w(nGrid-nResources))! allpcate size of vertical distribution of small pelagics     
      !allocate(self%id_smzoo_fish_c(nGrid-nResources), self%id_smzoo_fish_n(nGrid-nResources), self%id_smzoo_fish_p(nGrid-nResources))
      !allocate(self%id_lgzoo_fish_c(nGrid-nResources), self%id_lgzoo_fish_n(nGrid-nResources), self%id_lgzoo_fish_p(nGrid-nResources))
      allocate(self%id_excre_fish_n(nGrid-nResources), self%id_excre_fish_p(nGrid-nResources))
      allocate(self%id_respiration_fish_c(nGrid-nResources),self%id_respiration_fish_n(nGrid-nResources),self%id_respiration_fish_p(nGrid-nResources))
      allocate(self%id_feces_fish_c(nGrid-nResources), self%id_feces_fish_n(nGrid-nResources), self%id_feces_fish_p(nGrid-nResources))
      allocate(self%id_carcasses_fish_c(nGrid-nResources), self%id_carcasses_fish_n(nGrid-nResources), self%id_carcasses_fish_p(nGrid-nResources))

      allocate(self%id_fish_total_biomass(nGroups))
      


      
! theta_matrix and vertical_distribution_matrix preparation (MUST be done before creating child models)
      allocate (self%depth_array(self%n_depth))
      allocate (self%photic_depth_array(self%n_photic_depth))
      allocate (self%theta_matrix(nGrid,nGrid,self%n_depth,self%n_photic_depth)) ! 120 depth layers from 50 to 6000 m, interval 50 m
      allocate (self%vertical_distribution_matrix(int(self%max_depth)+1,nFGrid,self%n_depth,self%n_photic_depth))
      
      self%depth_array = [(real(i*self%depth_interval, rk), i=1, self%n_depth)] ! depth_array will contain: 50, 100, ..., 6000
      self%photic_depth_array = [(real(i*self%photic_depth_interval, rk), i=1, self%n_photic_depth)] ! photic depth array will contain: 1, 2, ..., 200
      self%theta_matrix=0._rk
      self%vertical_distribution_matrix=-999._rk  ! mark invalid data
         
      inquire(file='matrix.dat', exist=file_exists)

      if (file_exists) then

         print*, "Reading matrices from existing file..."
         open(unit=10, file='matrix.dat', form='unformatted', access='stream', status='old')
         read(10) dims
         if (dims(1)==size(self%theta_matrix,1) .and. dims(2)==size(self%theta_matrix,2) .and. dims(3)==size(self%theta_matrix,3) .and. dims(4)==size(self%theta_matrix,4)) then
            print *, "Dimensions match.", dims, ". Reading matrices..."
            read(10) self%theta_matrix
            read(10) self%vertical_distribution_matrix
            close(10)
            print *, "Matrices read successfully."
         else
            print *, "Dimension mismatch! Expected:", size(self%theta_matrix,1),size(self%theta_matrix,2),size(self%theta_matrix,3),size(self%theta_matrix,4), " Got:", dims, "."
            file_exists = .FALSE.
         end if      
      end if 
      
      if ( .NOT. file_exists ) then
         print*, "Preparing the theta matrix data file..."

         do i = 1, size(self%theta_matrix,3)
            depth = real(i*self%depth_interval,rk)
            do j = 1, size(self%theta_matrix,4)
               photic = real(j*self%photic_depth_interval,rk)
               call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,nStages,Tp,Tm,Tb,depth,photic,etamature,self%shelfdepth,visual,Fmax,etaF)
               self%theta_matrix(:,:,i,j)=theta
               
               ! Save vertical distribution (averaged day/night) - only fish (idxF:nGrid)
               self%vertical_distribution_matrix(1:size(depthDay,1),:,i,j) = (depthDay(:,idxF:nGrid) + depthNight(:,idxF:nGrid))/2._rk
               
               ! shallow water column on shelf: no mesopelagics and mid-water predators
               if(depth.le.self%shelfdepth) self%theta_matrix(ixStart(2):ixEnd(2),:,i,j)=0._rk
               if(depth.le.self%shelfdepth) self%theta_matrix(ixStart(4):ixEnd(4),:,i,j)=0._rk
            end do
         end do
         mort0 = mort0/365._rk/86400._rk
         mortF = mortF/365._rk/86400._rk

         print *, "Saving matrices to file..."
         dims(1:4)=[size(self%theta_matrix,1),size(self%theta_matrix,2),size(self%theta_matrix,3),size(self%theta_matrix,4)]
         open(unit=10, file='matrix.dat', form='unformatted', access='stream', status='replace')
         write(10) dims
         write(10) self%theta_matrix
         write(10) self%vertical_distribution_matrix
         close(10)
         print *, "Matrices saved successfully."
      end if

      print *, "Matrices (theta, vertical_distribution) ready for computation."





      
      do i = 1, nGroups
         write (i_str,'(i0)') i
         ! read functional type longnames from yaml
         call self%get_parameter(fft_long_name,  'fft_long_name'//trim(i_str),  units='', long_name='fish functional type long name',  default="fish_functional_type_"//trim(i_str))
         ! read functional type shortnames from yaml
         call self%get_parameter(fft_short_name, 'fft_short_name_'//trim(i_str), units='', long_name='fish functional type short name', default="fft_"//trim(i_str))
         ! read the initial value of each functional type from yaml
         call self%get_parameter(self%b_ini, 'b_ini_'//trim(i_str), 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  

            do j = ixStart(i)-nResources, ixEnd(i)-nResources
               write (j_str,'(i0)') j
               write (size_number_str,'(i0)') j+nResources-ixStart(i)+1
               ! assign the initial biomass to each fish size class of a functional type.
               !!change names!!
                call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(size_number_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(size_number_str)//'_biomass', initial_value=self%b_ini, minimum=0.0_rk)

               !if(i .eq. 2 .OR.i .eq. 4)then
               !   if(depth.le.shelfdepth) then 
               !      call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(size_number_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(size_number_str)//'_biomass', initial_value=0.0_rk, minimum=0.0_rk) 
               !   else
               !      call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(size_number_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(size_number_str)//'_biomass', initial_value=b_ini, minimum=0.0_rk)
               !   end if
               !else
               !   call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(size_number_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(size_number_str)//'_biomass', initial_value=b_ini, minimum=0.0_rk)   
               !end if
            end do
            
          call self%register_diagnostic_variable(self%id_fish_total_biomass(i), trim(fft_short_name)//'_totB', 'g m-2', trim(fft_short_name)//'_total_biomass')   
            
      end do 
      
      do i = 1, nGrid-nResources
         write (i_str,'(i0)') i
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_fish(i), scale_factor = gww_mmolC)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fish(i), scale_factor = gww_mmolN)
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_fish(i), scale_factor = gww_mmolP)
         
         ! Register vertical distribution for this fish size class
         call self%register_vertical_distribution(self%id_fish_w(i),'fish_'//trim(i_str))
         
         ! Create FEISTY-controlled vertical distribution model
         allocate(feisty_vertical_distribution)
         call self%add_child(feisty_vertical_distribution, 'feisty_vertical_distribution_'//trim(i_str), configunit=-1)
         call self%request_coupling(self%id_fish_w(i), 'feisty_vertical_distribution_'//trim(i_str)//'/'//'w')
         
         ! Set pointers to parent FEISTY model's data
         feisty_vertical_distribution%vertical_distribution_matrix => self%vertical_distribution_matrix
         feisty_vertical_distribution%depth_array => self%depth_array
         feisty_vertical_distribution%photic_depth_array => self%photic_depth_array
         feisty_vertical_distribution%fish_index = i
         
         !call self%register_state_dependency(self%id_smzoo_fish_c(i), 'smzoo_c', 'mmol C m-2', 'depth-integrated small zooplankton carbon')
         !call self%register_state_dependency(self%id_smzoo_fish_n(i), 'smzoo_n', 'mmol N m-2', 'depth-integrated small zooplankton nitrogen')
         !call self%register_state_dependency(self%id_smzoo_fish_p(i), 'smzoo_p', 'mmol P m-2', 'depth-integrated small zooplankton phosphorus')   
         !
         !call self%register_state_dependency(self%id_lgzoo_fish_c(i), 'lgzoo_c', 'mmol C m-2', 'depth-integrated large zooplankton carbon')
         !call self%register_state_dependency(self%id_lgzoo_fish_n(i), 'lgzoo_n', 'mmol N m-2', 'depth-integrated large zooplankton nitrogen')
         !call self%register_state_dependency(self%id_lgzoo_fish_p(i), 'lgzoo_p', 'mmol P m-2', 'depth-integrated large zooplankton phosphorus')   
         
         !call self%register_state_dependency(self%id_excre_fish_n(i),'excretion_n', 'mmol N m-2', 'excretion nitrogen')
         !call self%register_state_dependency(self%id_excre_fish_p(i),'excretion_p', 'mmol C m-2', 'excretion phosphorus')
         
         call self%register_state_dependency(self%id_respiration_fish_c(i),'respiration_c', 'mmol C m-2', 'respiration carbon')
         call self%register_state_dependency(self%id_respiration_fish_n(i),'respiration_n', 'mmol N m-2', 'respiration nitrogen')
         call self%register_state_dependency(self%id_respiration_fish_p(i),'respiration_p', 'mmol P m-2', 'respiration phosphorus')
         
         call self%register_state_dependency(self%id_feces_fish_c(i),'feces_c', 'mmol C m-2', 'feces carbon')
         call self%register_state_dependency(self%id_feces_fish_n(i),'feces_n', 'mmol N m-2', 'feces nitrogen')
         call self%register_state_dependency(self%id_feces_fish_p(i),'feces_p', 'mmol P m-2', 'feces phosphorus')
         
         call self%register_state_dependency(self%id_carcasses_fish_c(i),'carcasses_c', 'mmol C m-2', 'carcasses carbon')
         call self%register_state_dependency(self%id_carcasses_fish_n(i),'carcasses_n', 'mmol N m-2', 'carcasses nitrogen')
         call self%register_state_dependency(self%id_carcasses_fish_p(i),'carcasses_p', 'mmol P m-2', 'carcasses phosphorus')
         
         !call self%request_mapped_coupling_to_model(self%id_smzoo_fish_c(i), 'small_zooplankton_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         !call self%request_mapped_coupling_to_model(self%id_smzoo_fish_n(i), 'small_zooplankton_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         !call self%request_mapped_coupling_to_model(self%id_smzoo_fish_p(i), 'small_zooplankton_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
         !call self%couplings%set_string('small_zooplankton_fish_'//trim(i_str), "small_zooplankton")
         !
         !call self%request_mapped_coupling_to_model(self%id_lgzoo_fish_c(i), 'large_zooplankton_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         !call self%request_mapped_coupling_to_model(self%id_lgzoo_fish_n(i), 'large_zooplankton_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         !call self%request_mapped_coupling_to_model(self%id_lgzoo_fish_p(i), 'large_zooplankton_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
         !call self%couplings%set_string('large_zooplankton_fish_'//trim(i_str), "large_zooplankton")
       
        !call self%request_mapped_coupling_to_model(self%id_excre_fish_n(i), 'excretion_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         !call self%request_mapped_coupling_to_model(self%id_excre_fish_p(i), 'excretion_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))  
        ! call self%couplings%set('excretion_fish_'//trim(i_str), "excretion")
         
         call self%request_mapped_coupling_to_model(self%id_respiration_fish_c(i),'respiration_fish_'//trim(i_str),standard_variables%total_carbon, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_respiration_fish_n(i),'respiration_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         call self%request_mapped_coupling_to_model(self%id_respiration_fish_p(i),'respiration_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
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
      call self%register_state_dependency(self%id_det, 'detritus', 'mmol N m-3', 'detritus reaching the bottom for driving benthos')
      call self%register_state_dependency(self%id_nut, 'nutrient', 'mmol N m-3', 'nutrient in the bottom layer')
      
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_benthos, scale_factor = gww_mmolC)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_benthos, scale_factor = gww_mmolN)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_benthos, scale_factor = gww_mmolP)
      
       ! Depth-integrated dependencies
       call self%register_dependency(self%id_smzoo_c, 'small_zooplankton_c', 'mmol C m-2', 'depth-integrated small zooplankton carbon')
       call self%register_dependency(self%id_smzoo_n, 'small_zooplankton_n', 'mmol N m-2', 'depth-integrated small zooplankton nitrogen')
       call self%register_dependency(self%id_smzoo_p, 'small_zooplankton_p', 'mmol P m-2', 'depth-integrated small zooplankton phosphorus') 
       call self%register_dependency(self%id_lgzoo_c, 'large_zooplankton_c', 'mmol C m-2', 'depth-integrated large zooplankton carbon')
       call self%register_dependency(self%id_lgzoo_n, 'large_zooplankton_n', 'mmol N m-2', 'depth-integrated large zooplankton nitrogen')
       call self%register_dependency(self%id_lgzoo_p, 'large_zooplankton_p', 'mmol P m-2', 'depth-integrated large zooplankton phosphorus') 
       
       call self%register_vertical_distribution(self%dummy_w,'dummy')
       allocate(depth_distribution)
       call self%add_child(depth_distribution, 'dummy')
       call self%request_coupling(self%dummy_w, 'dummy'//'/'//'w')
       
       call self%request_mapped_coupling_to_model(self%id_smzoo_c, 'small_zooplankton', standard_variables%total_carbon, id_w=self%dummy_w)!, average=.true.)
       call self%request_mapped_coupling_to_model(self%id_smzoo_n, 'small_zooplankton', standard_variables%total_nitrogen, id_w=self%dummy_w)!, average=.true.)
       call self%request_mapped_coupling_to_model(self%id_smzoo_p, 'small_zooplankton', standard_variables%total_phosphorus, id_w=self%dummy_w)!, average=.true.)
       !call self%couplings%set_string('small_zooplankton', "small_zooplankton")
       call self%request_mapped_coupling_to_model(self%id_lgzoo_c, 'large_zooplankton', standard_variables%total_carbon, id_w=self%dummy_w)!, average=.true.)
       call self%request_mapped_coupling_to_model(self%id_lgzoo_n, 'large_zooplankton', standard_variables%total_nitrogen, id_w=self%dummy_w)!, average=.true.)
       call self%request_mapped_coupling_to_model(self%id_lgzoo_p, 'large_zooplankton', standard_variables%total_phosphorus, id_w=self%dummy_w)!, average=.true.)
       !call self%couplings%set_string('large_zooplankton', "large_zooplankton")
       
       !???????????????????????????????????????????????????????????????????????
       !???????????????????????????????????????????????????????????????????????
       !???????????????????????????????????????????????????????????????????????
       call self%register_dependency(self%id_kpar, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
       !======================================================================
       
       
       call self%register_dependency(self%id_temp, standard_variables%temperature)
       call self%register_dependency(self%id_temp_vertmean_100m, vertical_mean(self%id_temp, maximum_depth=100._rk))
       call self%register_dependency(self%id_temp_vertmean_500to1500m, vertical_mean(self%id_temp, minimum_depth=500._rk, maximum_depth=1500._rk))
       call self%register_diagnostic_variable(self%id_temp_vertmean_100m_diag, 'temp_vertmean_100m', 'degree_C', 'vertical mean temperature above 100 m')
       call self%register_dependency(self%id_bottom_depth, standard_variables%bottom_depth)
       
       call self%register_mapped_model_dependency(self%id_smzoo_int, 'small_zooplankton', proportional_change=.true., domain=domain_bottom,id_w=self%dummy_w) ! check default average
       call self%register_mapped_model_dependency(self%id_lgzoo_int, 'large_zooplankton', proportional_change=.true., domain=domain_bottom,id_w=self%dummy_w)
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
      real(rk) :: ingestion_c, ingestion_n, ingestion_p, prey_loss_rate, p1,p2, net_growth
      integer  :: iGroup, i, istate, closest_bottom_depth_idx,closest_photic_depth_idx
      real(rk),dimension(ixEnd(1)-ixStart(1)+1)  :: smpel
      real(rk),dimension(ixEnd(2)-ixStart(2)+1)  :: lgpel
      real(rk),dimension(ixEnd(3)-ixStart(3)+1)  :: dem
      real(rk),dimension(nGrid-nResources)       :: fish
      real(rk)                                   :: smzoo_c, smzoo_n, smzoo_p, lgzoo_c, lgzoo_n, lgzoo_p
      real(rk)                                   :: excessive_c,excessive_n,excessive_p
      real(rk)                                   :: excessive_c_smzoo,excessive_n_smzoo, excessive_p_smzoo
      real(rk)                                   :: excessive_c_lgzoo,excessive_n_lgzoo, excessive_p_lgzoo
      real(rk)                                   :: zoo1, zoo2, benthos1, benthos2
      real(rk)                                   :: det_bot, det_bot_flux_gww, benthos_g_gww, benthos_loss_gww
      real(rk),dimension(nGrid)                  :: uin, dudt, mortpred_contri_zoo1, mortpred_contri_zoo2
      real(rk)                                   :: temp_100m,temp_bottom, temp_500to1500m, bottom_depth
      real(rk)                                   :: kpar,photic_depth
      real(rk), allocatable, dimension(:,:)      :: vertical_dist  ! averaged vertical distribution (day+night)/2 (xrange, nFGrid - fish only)
      
      _BOTTOM_LOOP_BEGIN_

         ! Get depth-integrated predator biomass       
         do i = 1, nGrid-nResources
         _GET_BOTTOM_(self%id_fish(i), fish(i))
         end do      
         !print*,(sum(fish))
         
         ! get depth-integrated small zooplankton (mmol N m-2)
         _GET_BOTTOM_(self%id_smzoo_c, smzoo_c) 
         _GET_BOTTOM_(self%id_smzoo_n, smzoo_n)
         _GET_BOTTOM_(self%id_smzoo_p, smzoo_p)
         ! get depth-integrated large zooplankton (mmol N m-2)
         _GET_BOTTOM_(self%id_lgzoo_c, lgzoo_c)
         _GET_BOTTOM_(self%id_lgzoo_n, lgzoo_n)
         _GET_BOTTOM_(self%id_lgzoo_p, lgzoo_p)
         !print*,smzoo_c/smzoo_n
         !print*,smzoo_c/smzoo_p
         
         zoo1=0._rk
         zoo2=0._rk
         ! mmol N/m2 to gww/m2
         zoo1 = minval([smzoo_c/gww_mmolC, smzoo_n/gww_mmolN, smzoo_p/gww_mmolP], mask=[smzoo_c/gww_mmolC .ne. 0._rk, smzoo_n/gww_mmolN .ne. 0._rk, smzoo_p/gww_mmolP .ne. 0._rk])
         zoo2 = minval([lgzoo_c/gww_mmolC, lgzoo_n/gww_mmolN, lgzoo_p/gww_mmolP], mask=[lgzoo_c/gww_mmolC .ne. 0._rk, lgzoo_n/gww_mmolN .ne. 0._rk, lgzoo_p/gww_mmolP .ne. 0._rk])
         !print*,smzoo_c/gww_mmolC,smzoo_n/gww_mmolN,smzoo_p/gww_mmolP
         
         ! split two size classes
         zoo1 = zoo1/2._rk
         zoo2 = zoo2/2._rk
         
         !get benthos (gww m-2)
         _GET_BOTTOM_(self%id_benthos,benthos1)
         _GET_BOTTOM_(self%id_benthos,benthos2)
         !print*,benthos1
         
         ! Depth-averaged environmental dependencies
         _GET_BOTTOM_(self%id_bottom_depth, bottom_depth) 
         _GET_BOTTOM_(self%id_temp_vertmean_100m, temp_100m)
         _GET_(self%id_temp, temp_bottom)
         _GET_BOTTOM_(self%id_temp_vertmean_500to1500m,temp_500to1500m)
         if(bottom_depth .lt. 500._rk) temp_500to1500m = temp_bottom
         _GET_(self%id_kpar, kpar)
         if (kpar > 1.0e-10_rk) then
            photic_depth = log(0.01_rk)/(-kpar) ! get the photic zone depth from the light attenuation coeffcient
         else
            ! If kpar is zero or near-zero (no attenuation), entire column is photic
            photic_depth = bottom_depth
         end if
         !print*, 'photic depth', photic_depth
         !print*,temp_100m
         _SET_BOTTOM_DIAGNOSTIC_(self%id_temp_vertmean_100m_diag, temp_100m)
         

         !call updateTemp(temp_100m, temp_bottom, bottom_depth, [1,2],2,[3],1)
         !call set2vec
         !mort0 = mort0/365._rk/86400._rk
         !mortF = mortF/365._rk/86400._rk
         !call setupbasic(100._rk, 100._rk, 100._rk, -1._rk, bottom_depth, temp_100m, temp_bottom)

         !update theta and vertical distribution
         closest_bottom_depth_idx = minloc(abs(self%depth_array - bottom_depth), dim=1)
         closest_photic_depth_idx= minloc(abs(self%photic_depth_array - photic_depth), dim=1)
         !print*, closest_depth_idx
         theta= self%theta_matrix(:,:,closest_bottom_depth_idx,closest_photic_depth_idx)
         
         ! Retrieve vertical distribution from matrix (only fish)
         if (allocated(vertical_dist)) deallocate(vertical_dist)
         allocate(vertical_dist(int(bottom_depth)+1, nFGrid))
         vertical_dist = self%vertical_distribution_matrix(1:int(bottom_depth)+1,:,closest_bottom_depth_idx,closest_photic_depth_idx)
         
         dvm = photic_depth + 500._dp ! 650._dp
         if (bottom_depth .lt. (photic_depth + 500._dp)) dvm = bottom_depth 
         if (bottom_depth .le. self%shelfdepth) dvm = 0._dp

         !update temperature effect
         call updateTempV2(temp_100m, temp_500to1500m, temp_bottom, dvm, bottom_depth, photic_depth, ixmedium, ixlarge)
         do iGroup = 1, nGroups
             group(iGroup)%spec%V=group(iGroup)%spec%Vsave*fTempV(ixStart(iGroup):ixEnd(iGroup))
             group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmaxsave*fTempV(ixStart(iGroup):ixEnd(iGroup))
             group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolismsave*fTempmV(ixStart(iGroup):ixEnd(iGroup))
         end do
         call set2vec
         mort0 = mort0/365._rk/86400._rk
         mortF = mortF/365._rk/86400._rk

         
         uin= [zoo1,zoo2 ,benthos1,0._rk ,fish]
         !print*,uin
         !uin= [100._rk,100._rk,5._rk,0._rk,fish]
         call calcderivatives(uin, dudt)
         !print*,dudt
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
         
         !small zooplankton
         do istate = 1, size(self%id_smzoo_int%bottom_state)
            _GET_BOTTOM_(self%id_smzoo_int%bottom_state(istate), p1)
            _ADD_BOTTOM_SOURCE_(self%id_smzoo_int%bottom_state(istate), dudt(1)/zoo1 * p1/2._rk)
         end do     
         !large zooplankton
         do istate = 1, size(self%id_lgzoo_int%bottom_state)
            _GET_BOTTOM_(self%id_lgzoo_int%bottom_state(istate), p2)
            _ADD_BOTTOM_SOURCE_(self%id_lgzoo_int%bottom_state(istate), dudt(2)/zoo2 * p2/2._rk)
         end do 
            !print*,dudt(1)/zoo1 * p1/2._rk + dudt(2)/zoo2 * p2/2._rk
            !print*,dudt(1)* gww_mmolN+dudt(2)* gww_mmolN
                
         !temperatory put here
         excessive_c_smzoo = smzoo_c/gww_mmolC/2._rk * dudt(1)/zoo1 - dudt(1) ! gww biomass * specific loss rate - loss
         excessive_n_smzoo = smzoo_n/gww_mmolN/2._rk * dudt(1)/zoo1 - dudt(1) 
         excessive_p_smzoo = smzoo_p/gww_mmolP/2._rk * dudt(1)/zoo1 - dudt(1)
         excessive_c = MAX(excessive_c_smzoo,0._rk)
         excessive_n = MAX(excessive_n_smzoo,0._rk)
         excessive_p = MAX(excessive_p_smzoo,0._rk)
         excessive_c_lgzoo = lgzoo_c/gww_mmolC/2._rk * dudt(2)/zoo2 - dudt(2)
         excessive_n_lgzoo = lgzoo_n/gww_mmolN/2._rk * dudt(2)/zoo2 - dudt(2)
         excessive_p_lgzoo = lgzoo_p/gww_mmolP/2._rk * dudt(2)/zoo2 - dudt(2)
         excessive_c = excessive_c + MAX(excessive_c_lgzoo,0._rk)
         excessive_n = excessive_c + MAX(excessive_n_lgzoo,0._rk)
         excessive_p = excessive_c + MAX(excessive_p_lgzoo,0._rk)
         
         
         
         do i = 1, nGrid-nResources
         ! fish dynamics
         _ADD_BOTTOM_SOURCE_(self%id_fish(i), dudt(i+nResources) )!gww/m2
         ! zooplankton dynamics
         !_ADD_BOTTOM_SOURCE_(self%id_smzoo_fish_c(i), -mortpred_contri_zoo1(i+nResources) * gww_mmolC)!gww/m2 to mmol C/m2
         !_ADD_BOTTOM_SOURCE_(self%id_smzoo_fish_n(i), -mortpred_contri_zoo1(i+nResources) * gww_mmolN)!/0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol N/m2
         !_ADD_BOTTOM_SOURCE_(self%id_smzoo_fish_p(i), -mortpred_contri_zoo1(i+nResources) * gww_mmolP)!gww/m2 to mmol P/m2
         !
         !_ADD_BOTTOM_SOURCE_(self%id_lgzoo_fish_c(i), -mortpred_contri_zoo2(i+nResources) * gww_mmolC)
         !_ADD_BOTTOM_SOURCE_(self%id_lgzoo_fish_n(i), -mortpred_contri_zoo2(i+nResources) * gww_mmolN)!/0.01201_rk/9_rk *16._rk/106._rk)
         !_ADD_BOTTOM_SOURCE_(self%id_lgzoo_fish_p(i), -mortpred_contri_zoo2(i+nResources) * gww_mmolP)
         
         ! excretion
         !_ADD_BOTTOM_SOURCE_(self%id_excre_fish_n(i), respiration(i+nResources) * gww_mmolN)!* CN /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol N/m2
         !_ADD_BOTTOM_SOURCE_(self%id_excre_fish_p(i), respiration(i+nResources) * gww_mmolP)!gww/m2 to mmol P/m2
         ! respiration
         _ADD_BOTTOM_SOURCE_(self%id_respiration_fish_c(i), respiration(i+nResources) * gww_mmolC)!* gwwC /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol C/m2
         _ADD_BOTTOM_SOURCE_(self%id_respiration_fish_n(i), respiration(i+nResources) * gww_mmolN)
         _ADD_BOTTOM_SOURCE_(self%id_respiration_fish_p(i), respiration(i+nResources) * gww_mmolP)
         ! feces
         _ADD_BOTTOM_SOURCE_(self%id_feces_fish_c(i), feces(i+nResources) * gww_mmolC)!gww/m2 to mmol C/m2
         _ADD_BOTTOM_SOURCE_(self%id_feces_fish_n(i), feces(i+nResources) * gww_mmolN)!* CN /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol N/m2
         _ADD_BOTTOM_SOURCE_(self%id_feces_fish_p(i), feces(i+nResources) * gww_mmolP)!gww/m2 to mmol P/m2
         ! carcasses
         _ADD_BOTTOM_SOURCE_(self%id_carcasses_fish_c(i), carcasses(i+nResources) * gww_mmolC)!gww/m2 to mmol C/m2
         _ADD_BOTTOM_SOURCE_(self%id_carcasses_fish_n(i), carcasses(i+nResources) * gww_mmolN)!* CN /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol N/m2
         _ADD_BOTTOM_SOURCE_(self%id_carcasses_fish_p(i), carcasses(i+nResources) * gww_mmolP)!gww/m2 to mmol P/m2
         end do
         
         ! benthos dynamics  
         _GET_(self%id_det,det_bot)! mmol N m-3
         det_bot_flux_gww=det_bot * 5.0_rk/86400.0_rk /gww_mmolN !det_bot*w_d   mmol N m-3 to mmol N m-3 * m s-1 to gww/m2 s-1
         i = 3!,3!, 4
         benthos_g_gww = 0.1_rk*det_bot_flux_gww*(1-benthos1/80)!gww/m2/s
         _ADD_BOTTOM_SOURCE_(self%id_benthos, benthos_g_gww + dudt(i) )!gww/m2  dRdt(3) = rr(3)*(1-R(3)/K(3)) - mortRes(3)*R(3) K is 80
         _ADD_BOTTOM_FLUX_(self%id_det, -det_bot_flux_gww*gww_mmolN) ! gww/m2 s-1 to mmol N/m2 s-1

         benthos_loss_gww = 0.1_rk*det_bot_flux_gww*(1-(1-benthos1/80)) !gww, convert below
         !print*, det_bot_flux_gww
         !print*, (0.9_rk *det_bot_flux_gww + benthos_loss_gww) +benthos_g_gww
         _ADD_BOTTOM_FLUX_(self%id_nut, (0.9_rk *det_bot_flux_gww + benthos_loss_gww) * gww_mmolN) !
         
         
         do i=1,nGroups            
            _SET_BOTTOM_DIAGNOSTIC_(self%id_fish_total_biomass(i), totBiomass(i))
         end do
         
         
         !end do
         
         !mass conservation check
         !fish
         !print*,-sum(mortpred_contri_zoo1(5:nGrid)) * gww_mmolN-sum(mortpred_contri_zoo2(5:nGrid)) * gww_mmolN + dudt(3)* gww_mmolN+ &
         !   & sum(dudt(5:nGrid))* gww_mmolN + sum(excretion + respiration +carcasses + feces)* gww_mmolN 
         !print*,sum(dudt(1:nGrid))+ sum(excretion + respiration +carcasses + feces)
         !print*, sum(dudt(1:4))* gww_mmolN+sum(dudt(5:nGrid))* gww_mmolN +sum(excretion + respiration +carcasses + feces)* gww_mmolN
          !print*, dudt(1)/zoo1 * p1/2._rk +dudt(2)/zoo2 * p2/2._rk + sum(dudt(3:4))* gww_mmolN+sum(dudt(5:nGrid))* gww_mmolN +sum(excretion + respiration +carcasses + feces)* gww_mmolN
         
         !benthos
         !print*, 0.1_rk*det_bot_flux_gww - benthos_g_gww - benthos_loss_gww
         !print*, -det_bot_flux_gww+(0.9_rk *det_bot_flux_gww + benthos_loss_gww)+benthos_g_gww

         
         ! Calculate ingested fluxes of different chemical elements
         ! Predator population growth will be based on the most limiting of these
         !ingestion_c = self%clearance_rate * c * prey_c
         !ingestion_n = self%clearance_rate * c * prey_n
         !ingestion_p = self%clearance_rate * c * prey_p
         !net_growth = min(ingestion_c, ingestion_n / CN, ingestion_p / CP) - self%mortality * c

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
         !_ADD_BOTTOM_SOURCE_(self%id_waste_n, ingestion_n - net_growth * CN)
         !_ADD_BOTTOM_SOURCE_(self%id_waste_p, ingestion_p - net_growth * CP)

         ! Save diagnostics
         !_SET_BOTTOM_DIAGNOSTIC_(self%id_net_growth, net_growth * 86400.0_rk)
         !_SET_BOTTOM_DIAGNOSTIC_(self%id_prey_loss_rate, prey_loss_rate * 86400.0_rk)

         
      _BOTTOM_LOOP_END_
   end subroutine   
   
   ! Add model subroutines here.

end module
