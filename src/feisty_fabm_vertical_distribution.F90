#include "fabm_driver.h"

module feisty_fabm_vertical_distribution
   ! This module provides FEISTY-calculated vertical distribution weights to FABM
   ! It reads from FEISTY's pre-calculated vertical_distribution_matrix_centers each timestep
   
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
      
      ! Pointers to parent FEISTY model's data (set during initialization)
      real(rk), pointer :: vertical_distribution_matrix_centers(:,:,:) => null()  ! (depth_center, nFGrid, depth) - contains center values
      real(rk), pointer :: depth_array(:) => null()
      integer  :: fish_index  ! Which fish size class this distribution is for
      
   contains
      procedure :: initialize => feisty_vertical_distribution_initialize
      procedure :: do_column  => feisty_vertical_distribution_do_column
   end type

contains

   subroutine feisty_vertical_distribution_initialize(self, configunit)
      class (type_feisty_vertical_distribution), intent(inout), target :: self
      integer,                   intent(in)            :: configunit
      
      call self%register_implemented_routines((/source_do_column/))
      
      ! Get which fish index this distribution is for
      call self%get_parameter(self%fish_index, 'fish_index', '', 'fish size class index', default=1)
      
      ! Register diagnostic variable
      call self%register_diagnostic_variable(self%id_w, 'w', '1', 'vertical distribution weights', source=source_do_column)
      
      ! Register dependencies
      call self%register_dependency(self%id_z, standard_variables%depth)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)
      call self%register_dependency(self%id_bottom_depth, standard_variables%bottom_depth)
   end subroutine feisty_vertical_distribution_initialize
   
   
   subroutine feisty_vertical_distribution_do_column(self, _ARGUMENTS_DO_COLUMN_)
      ! This routine is called once per water column
      ! It assigns weights to layers based on pre-calculated vertical distribution
      
      class (type_feisty_vertical_distribution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_
      
      real(rk) :: z_center, h, z_top, z_bottom, bottom_depth
      real(rk) :: weight_sum, weight
      real(rk) :: max_depth
      real(rk) :: interval_lower, interval_upper, coverage_fraction
      real(rk) :: scale_factor, profile_depth
      real(rk) :: z_top_feisty, z_bottom_feisty
      real(rk) :: total_assigned_weight  ! Sum of all assigned weights for testing
      integer  :: depth_idx
      integer  :: closest_bottom_depth_idx
      integer  :: z_start_idx, z_end_idx
      
      ! Check if matrix data is available
      ! if (.not. associated(self%vertical_distribution_matrix_centers)) then
      !    ! No data available - use uniform distribution
      !    _DOWNWARD_LOOP_BEGIN_
      !       _SET_DIAGNOSTIC_(self%id_w, 1.0_rk)
      !    _DOWNWARD_LOOP_END_
      !    return
      ! end if
      
      ! Get bottom depth (same for all layers in this column)
      _GET_HORIZONTAL_(self%id_bottom_depth, bottom_depth)
      ! Find the closest pre-computed depth profile
      closest_bottom_depth_idx = minloc(abs(self%depth_array - bottom_depth), dim=1)
      ! Get the FEISTY profile depth for this index
      profile_depth = self%depth_array(closest_bottom_depth_idx)
      ! Calculate scale factor to compress FEISTY distribution to FABM water column
      ! This ensures all FEISTY weights (0 to profile_depth) are distributed into FABM (0 to bottom_depth)
      scale_factor = bottom_depth / profile_depth
      ! Use full FEISTY profile depth for calculations
      ! Matrix contains center values: index 1 = 0.5m, index 2 = 1.5m, ..., index N = (N-0.5)m
      max_depth = min(profile_depth, real(size(self%vertical_distribution_matrix_centers, 1), rk))
      
      ! Initialize accumulator for testing weight summation
      total_assigned_weight = 0.0_rk
      
      ! Calculate vertical distribution weights for each layer
      _DOWNWARD_LOOP_BEGIN_
         ! Get layer properties from FABM
         _GET_(self%id_z, z_center)  ! Layer center depth (m)
         _GET_(self%id_h, h)          ! Layer thickness (m)
         
         ! Calculate layer boundaries (in meters) - FABM coordinates
         z_top = z_center - 0.5_rk * h
         z_bottom = z_center + 0.5_rk * h
         
         ! Scale FABM depths to FEISTY coordinates for mapping
         ! This compresses the FEISTY distribution (0 to profile_depth) into FABM water column (0 to bottom_depth)
         z_top_feisty = z_top / scale_factor
         z_bottom_feisty = z_bottom / scale_factor
         
         ! Universal approach: calculate weight by summing contributions from each meter interval
         ! Matrix index N corresponds to center at (N - 0.5) meters, representing interval [N-1, N]
         ! Find which meter intervals overlap with this layer in FEISTY coordinates
         ! Use floor() to find the interval containing each boundary: interval index = floor(z) + 1
         z_start_idx = max(1, int(floor(z_top_feisty)) + 1)  ! First interval index (interval containing z_top_feisty)
         z_end_idx = min(int(floor(max_depth)) + 1, int(floor(z_bottom_feisty)) + 1)  ! Last interval index (interval containing z_bottom_feisty)
         
         if (z_end_idx >= z_start_idx .and. z_start_idx >= 1) then
            weight_sum = 0.0_rk
            
            ! Loop through each meter interval that overlaps with this layer
            do depth_idx = z_start_idx, z_end_idx
               if (depth_idx >= 1 .and. depth_idx <= size(self%vertical_distribution_matrix_centers, 1)) then
                  ! Calculate coverage fraction for this interval [depth_idx-1, depth_idx]
                  interval_lower = real(depth_idx - 1, rk)  ! Lower boundary (e.g., 0m for idx 1)
                  interval_upper = real(depth_idx, rk)      ! Upper boundary (e.g., 1m for idx 1)
                  
                  ! Calculate how much of this interval is covered by the layer (in FEISTY coordinates)
                  ! For partial intervals (first and last), coverage < 1.0
                  ! For fully covered intervals (middle), coverage = 1.0
                  coverage_fraction = (min(z_bottom_feisty, interval_upper) - max(z_top_feisty, interval_lower)) / 1.0_rk
                  coverage_fraction = max(0.0_rk, min(1.0_rk, coverage_fraction))  ! Clamp to [0, 1]
                  
                  ! Add scaled weight: coverage_fraction * full_weight
                  weight_sum = weight_sum + &
                     coverage_fraction * self%vertical_distribution_matrix_centers(depth_idx, self%fish_index, closest_bottom_depth_idx)
               end if
            end do
            
            ! Integrated weight for this layer
            weight = weight_sum
         else
            ! Out of range
            weight = 0.0_rk
         end if

         ! Set the diagnostic weight for this layer
         _SET_DIAGNOSTIC_(self%id_w, weight)
         
         ! Accumulate weight for testing summation
         total_assigned_weight = total_assigned_weight + weight
      _DOWNWARD_LOOP_END_
      
      ! Test: Print the sum of all assigned weights (should be 1.0)
      !print*, 'Fish index:', self%fish_index
      !print*, 'Total assigned weight sum after scaling:', total_assigned_weight
      !print*, 'Original FEISTY weight sum (before scaling):', sum(self%vertical_distribution_matrix_centers(1:int(max_depth), self%fish_index, closest_bottom_depth_idx))
      !print*, 'Scale factor:', scale_factor, 'Bottom depth:', bottom_depth, 'Profile depth:', profile_depth
   end subroutine feisty_vertical_distribution_do_column

end module feisty_fabm_vertical_distribution

