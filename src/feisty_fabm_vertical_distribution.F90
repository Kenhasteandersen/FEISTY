#include "fabm_driver.h"

module feisty_fabm_vertical_distribution
   ! This module calculates FEISTY vertical distribution weights for the layers of the host model
   ! It maps FEISTY's pre-calculated 1 m profiles (vertical_distribution_matrix_centers) onto the host grid

   use fabm_types

   implicit none

   private

   public feisty_layer_weights

contains

   subroutine feisty_layer_weights(profiles, z_center, h, scale_factor, max_depth, w)
      ! Weights of one host layer for all fish size classes (see type_vertical_depth_range in depth_mapping.F90)
      real(rk), intent(in)  :: profiles(:,:)   ! (depth_center, fish size class); index N = centre at (N - 0.5) m
      real(rk), intent(in)  :: z_center, h     ! layer centre depth and thickness (m)
      real(rk), intent(in)  :: scale_factor    ! host bottom depth / FEISTY profile depth
      real(rk), intent(in)  :: max_depth       ! deepest FEISTY profile depth used (m)
      real(rk), intent(out) :: w(:)            ! weights for this layer, one per fish size class

      real(rk) :: z_top_feisty, z_bottom_feisty, coverage_fraction
      integer  :: depth_idx, z_start_idx, z_end_idx

      ! Scale host depths to FEISTY coordinates, compressing the FEISTY profile (0 to profile_depth) into the water column
      z_top_feisty    = (z_center - 0.5_rk * h) / scale_factor
      z_bottom_feisty = (z_center + 0.5_rk * h) / scale_factor

      ! Find the 1 m intervals that overlap with this layer (interval index = floor(z) + 1)
      z_start_idx = max(1, int(floor(z_top_feisty)) + 1)
      z_end_idx   = min(int(floor(max_depth)) + 1, int(floor(z_bottom_feisty)) + 1, size(profiles, 1))

      ! Sum the contributions of the overlapping intervals, weighted by the fraction of each interval covered by the layer
      w = 0.0_rk
      do depth_idx = z_start_idx, z_end_idx
         coverage_fraction = min(z_bottom_feisty, real(depth_idx, rk)) - max(z_top_feisty, real(depth_idx - 1, rk))
         coverage_fraction = max(0.0_rk, min(1.0_rk, coverage_fraction))  ! Clamp to [0, 1]
         w = w + coverage_fraction * profiles(depth_idx, :)
      end do

      ! Integrated weight for this layer
      w = w / h
   end subroutine feisty_layer_weights

end module feisty_fabm_vertical_distribution
