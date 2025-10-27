#include "fabm_driver.h"

module feisty_fabm_vertical_distribution
   ! This module provides FEISTY-calculated vertical distribution weights to FABM
   ! It reads from FEISTY's pre-calculated vertical_distribution_matrix each timestep
   
   use fabm_types
   
   implicit none
   
   private
   
   public type_feisty_vertical_distribution
   
   type, extends(type_base_model) :: type_feisty_vertical_distribution  ! See type_vertical_depth_range in depth_mapping.F90
      ! FABM diagnostic for weights
      type (type_diagnostic_variable_id) :: id_w   ! weights for this layer
      
      ! Dependencies for determining which distribution to use
      type (type_dependency_id)          :: id_z   ! depth of layer center (m)
      type (type_dependency_id)          :: id_h   ! layer thickness (m)
      type (type_horizontal_dependency_id) :: id_bottom_depth  ! bottom depth
      type (type_dependency_id)          :: id_kpar ! light attenuation coefficient
      
      ! Pointers to parent FEISTY model's data (set during initialization)
      real(rk), pointer :: vertical_distribution_matrix(:,:,:,:) => null()
      real(rk), pointer :: depth_array(:) => null()
      real(rk), pointer :: photic_depth_array(:) => null()
      integer  :: fish_index  ! Which fish size class this distribution is for
      
   contains
      procedure :: initialize => feisty_vertical_distribution_initialize
      procedure :: do         => feisty_vertical_distribution_do
      procedure :: do_column  => feisty_vertical_distribution_do_column
   end type

contains

   subroutine feisty_vertical_distribution_initialize(self, configunit)
      class (type_feisty_vertical_distribution), intent(inout), target :: self
      integer,                   intent(in)            :: configunit
      
      ! Register that this model has 'do' and 'do_column' routines
      ! do_column is preferred (optimal), but do() works as fallback for GOTM
      call self%register_implemented_routines((/source_do, source_do_column/))
      
      ! Get which fish index this distribution is for
      call self%get_parameter(self%fish_index, 'fish_index', '', 'fish size class index', default=1)
      
      ! Register diagnostic variable
      call self%register_diagnostic_variable(self%id_w, 'w', '1', 'vertical distribution weights')
      
      ! Register dependencies
      call self%register_dependency(self%id_z, standard_variables%depth)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)
      call self%register_dependency(self%id_bottom_depth, standard_variables%bottom_depth)
      call self%register_dependency(self%id_kpar, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
   end subroutine feisty_vertical_distribution_initialize
   
   
   subroutine feisty_vertical_distribution_do(self, _ARGUMENTS_DO_)
      ! Fallback routine for hosts that don't support do_column (e.g., GOTM)
      ! Uses per-layer kpar to estimate photic depth (approximation - not optimal)
      
      class (type_feisty_vertical_distribution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: z_center, h, z_top, z_bottom, bottom_depth, kpar, photic_depth
      real(rk) :: weight_sum, weight
      integer  :: z_start_m, z_end_m, depth_m, n_meters, max_depth_m
      integer  :: closest_bottom_depth_idx, closest_photic_depth_idx
      
      ! Check if matrix data is available
      if (.not. associated(self%vertical_distribution_matrix)) then
         ! No data available - use uniform distribution
         _LOOP_BEGIN_
            _SET_DIAGNOSTIC_(self%id_w, 1.0_rk)
         _LOOP_END_
         return
      end if
      
      ! Get bottom depth (same for all layers in this column)
      _GET_HORIZONTAL_(self%id_bottom_depth, bottom_depth)
      closest_bottom_depth_idx = minloc(abs(self%depth_array - bottom_depth), dim=1)
      max_depth_m = min(int(bottom_depth), size(self%vertical_distribution_matrix, 1) - 1)
      
      _LOOP_BEGIN_
         ! Get layer properties from FABM
         _GET_(self%id_z, z_center)  ! Layer center depth (m)
         _GET_(self%id_h, h)          ! Layer thickness (m)
         _GET_(self%id_kpar, kpar)    ! Light attenuation coefficient
         
         ! Estimate photic depth from this layer's kpar (approximation!)
         if (kpar > 1.0e-10_rk) then
            photic_depth = abs(log(0.01_rk)) / kpar  ! depth where I = 0.01 * I0
         else
            photic_depth = bottom_depth
         end if
         closest_photic_depth_idx = minloc(abs(self%photic_depth_array - photic_depth), dim=1)
         
         ! Calculate layer boundaries (in meters)
         z_top = z_center - 0.5_rk * h
         z_bottom = z_center + 0.5_rk * h
         
         ! Convert to integer indices for FEISTY's per-meter array (1-based from surface)
         z_start_m = max(1, int(floor(z_top)) + 1)
         z_end_m = min(max_depth_m, int(ceiling(z_bottom)))
         
         ! Calculate average weight over this layer from FEISTY's per-meter weights
         if (z_end_m >= z_start_m .and. z_start_m <= max_depth_m) then
            weight_sum = 0.0_rk
            n_meters = 0
            
            ! Sum FEISTY's per-meter weights over this FABM layer
            do depth_m = z_start_m, z_end_m
               if (depth_m >= 1 .and. depth_m <= size(self%vertical_distribution_matrix, 1)) then
                  weight_sum = weight_sum + &
                     self%vertical_distribution_matrix(depth_m, self%fish_index, &
                                                       closest_bottom_depth_idx, closest_photic_depth_idx)
                  n_meters = n_meters + 1
               end if
            end do
            
            ! Average weight for this layer
            if (n_meters > 0) then
               weight = weight_sum / real(n_meters, rk)
            else
               weight = 0.0_rk
            end if
         else
            ! Out of range
            weight = 0.0_rk
         end if
         
         ! Set the diagnostic weight for this layer
         _SET_DIAGNOSTIC_(self%id_w, weight)
      _LOOP_END_
   end subroutine feisty_vertical_distribution_do
   
   
   subroutine feisty_vertical_distribution_do_column(self, _ARGUMENTS_DO_COLUMN_)
      ! This routine is called once per water column
      ! It calculates photic depth by accumulating light extinction, then assigns weights to layers
      
      class (type_feisty_vertical_distribution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_
      
      real(rk) :: z_center, h, z_top, z_bottom, bottom_depth, kpar, photic_depth
      real(rk) :: weight_sum, weight, cumulative_extinction, depth_cumulative
      integer  :: z_start_m, z_end_m, depth_m, n_meters, max_depth_m
      integer  :: closest_bottom_depth_idx, closest_photic_depth_idx
      logical  :: photic_found
      
      ! Check if matrix data is available
      if (.not. associated(self%vertical_distribution_matrix)) then
         ! No data available - use uniform distribution
         _DOWNWARD_LOOP_BEGIN_
            _SET_DIAGNOSTIC_(self%id_w, 1.0_rk)
         _DOWNWARD_LOOP_END_
         return
      end if
      
      ! Get bottom depth (same for all layers in this column)
      _GET_HORIZONTAL_(self%id_bottom_depth, bottom_depth)
      closest_bottom_depth_idx = minloc(abs(self%depth_array - bottom_depth), dim=1)
      max_depth_m = min(int(bottom_depth), size(self%vertical_distribution_matrix, 1) - 1)
      
      ! First pass: Calculate photic depth by accumulating light extinction from surface to bottom
      ! Photic depth is where cumulative optical depth = abs(ln(0.01)) = abs(ln(I/I0))
      cumulative_extinction = 0.0_rk
      depth_cumulative = 0.0_rk
      photic_found = .false.
      photic_depth = bottom_depth  ! Default: entire column is photic
      
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_h, h)
         _GET_(self%id_kpar, kpar)
         
         ! Accumulate extinction through this layer: optical_depth += kpar * layer_thickness
         cumulative_extinction = cumulative_extinction + kpar * h
         depth_cumulative = depth_cumulative + h
         
         ! Check if we've reached 1% light level (optical depth = abs(ln(0.01)))
         if (.not. photic_found .and. cumulative_extinction >= abs(log(0.01_rk))) then
            ! Interpolate to find exact photic depth within this layer
            photic_depth = depth_cumulative - h + &
               (abs(log(0.01_rk)) - (cumulative_extinction - kpar * h)) / max(kpar, 1.0e-10_rk)
            photic_found = .true.
         end if
      _DOWNWARD_LOOP_END_
      
      ! Now we have the photic depth for this column
      closest_photic_depth_idx = minloc(abs(self%photic_depth_array - photic_depth), dim=1)
      
      ! Second pass: Calculate vertical distribution weights for each layer
      _DOWNWARD_LOOP_BEGIN_
         ! Get layer properties from FABM
         _GET_(self%id_z, z_center)  ! Layer center depth (m)
         _GET_(self%id_h, h)          ! Layer thickness (m)
         
         ! Calculate layer boundaries (in meters)
         z_top = z_center - 0.5_rk * h
         z_bottom = z_center + 0.5_rk * h
         
         ! Convert to integer indices for FEISTY's per-meter array (1-based from surface)
         z_start_m = max(1, int(floor(z_top)) + 1)
         z_end_m = min(max_depth_m, int(ceiling(z_bottom)))
         
         ! Calculate average weight over this layer from FEISTY's per-meter weights
         if (z_end_m >= z_start_m .and. z_start_m <= max_depth_m) then
            weight_sum = 0.0_rk
            n_meters = 0
            
            ! Sum FEISTY's per-meter weights over this FABM layer
            do depth_m = z_start_m, z_end_m
               if (depth_m >= 1 .and. depth_m <= size(self%vertical_distribution_matrix, 1)) then
                  weight_sum = weight_sum + &
                     self%vertical_distribution_matrix(depth_m, self%fish_index, &
                                                       closest_bottom_depth_idx, closest_photic_depth_idx)
                  n_meters = n_meters + 1
               end if
            end do
            
            ! Average weight for this layer
            if (n_meters > 0) then
               weight = weight_sum / real(n_meters, rk)
            else
               weight = 0.0_rk
            end if
         else
            ! Out of range
            weight = 0.0_rk
         end if
         
         ! Set the diagnostic weight for this layer
         _SET_DIAGNOSTIC_(self%id_w, weight)
      _DOWNWARD_LOOP_END_
   end subroutine feisty_vertical_distribution_do_column

end module feisty_fabm_vertical_distribution

