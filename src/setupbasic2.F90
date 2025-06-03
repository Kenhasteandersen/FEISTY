#include "fabm_driver.h"

module feisty_setupbasic2

   use globals
   use spectrum
   use fish
   use setup
   use fabm_types
   use fabm_particle
   use fabm_builtin_depth_mapping   
   
   implicit none

   private

   type, extends(type_depth_integrated_particle), public :: type_feisty_setupbasic2
      ! Add variable identifiers and parameters here.
      type (type_surface_state_variable_id),         allocatable :: id_u(:)
      
      
   contains
      procedure :: initialize
    !  procedure :: do_surface
      ! Reference model procedures here.
   end type type_feisty_setupbasic2

contains

   subroutine initialize(self, configunit)
      class (type_feisty_setupbasic2), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
      
      real(rk)           :: smz_ini, lgz_ini, smbent_ini, lgbent_ini, b_ini
      real(rk)           :: szprod, lzprod, bprodin, dfbot, depth, Tp, Tb
      real(rk)           :: etaMature, Fmax, etaF
      logical            :: bET
      integer :: i
      integer :: nStages, bETin
      character(len=10)  :: strindex
      
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
      ! for setupBasic2
      call self%get_parameter(nStages, 'nStages', '', 'size number of large fish functional groups', default=9)      
      call self%get_parameter(etaMature, 'etaMature', '-', 'coefficient determines the fish size with a 50% maturity level (etaMature*asymptotic size)', default=0.025_rk)    
      call self%get_parameter(Fmax, 'Fmax', 'yr-1', 'maximum fishing mortality', default=0._rk)  
      call self%get_parameter(etaF, 'etaF', '-', 'coefficient determines the fish size with a 50% fishing selectivity (etaF*asymptotic size)', default=0.05_rk)  
      call self%get_parameter(bET, 'bET', '', 'whether to enable effective temperature effects on the large demersal fish (1: on, 0: off)', default=.true.)             
      ! Register model parameters and variables here.
      
      !call setupbasic(szprod, lzprod, bprodin, dfbot, depth, Tp, Tb)
      if (bET == .true.) bETin = 1
      if (bET == .false.) bETin = 0
      call setupbasic2(szprod, lzprod, bprodin, dfbot, nStages, depth, Tp, Tb, etaMature,Fmax,etaF,bETin)
      
      ! Register state variables
      allocate(self%id_u(nGrid))
      call self%register_state_variable(self%id_u(1), 'u1', 'g m-2', 'biomass', initial_value=smz_ini, minimum=0.0_rk)
      call self%register_state_variable(self%id_u(2), 'u2', 'g m-2', 'biomass', initial_value=lgz_ini, minimum=0.0_rk)
      call self%register_state_variable(self%id_u(3), 'u3', 'g m-2', 'biomass', initial_value=smbent_ini, minimum=0.0_rk)
      call self%register_state_variable(self%id_u(4), 'u4', 'g m-2', 'biomass', initial_value=lgbent_ini, minimum=0.0_rk)      
      do i = idxF, nGrid
         write (strindex,'(i0)') i
         call self%register_state_variable(self%id_u(i), 'u'//trim(strindex), 'g m-2', 'biomass', initial_value=b_ini, minimum=0.0_rk)
               call self%set_variable_property(self%id_u(i), 'disable_transport', .true.)
      end do

      
   end subroutine initialize
   
   !subroutine do(self, _ARGUMENTS_DO_)
   !   class (type_feisty_setupbasic), intent(in) :: self
   !   _DECLARE_ARGUMENTS_DO_
   !
   !   real(rk) :: prey, clearance_rate
   !
   !   ! Enter spatial loops (if any).
   !   _LOOP_BEGIN_
   !
   !      ! Retrieve current (local) state variable values.
   !      _GET_(self%id_prey, prey)                        ! prey density
   !      if (self%use_external_clearance_rate) then
   !         _GET_(self%id_clearance_rate, clearance_rate) ! prescribed clearance rate
   !      else
   !         clearance_rate = self%clearance_rate
   !      end if
   !
   !      ! Set fluxes of pelagic variables.
   !      _ADD_SOURCE_(self%id_prey, -prey * clearance_rate)
   !      _ADD_SOURCE_(self%id_consumed_prey, prey * clearance_rate)
   !
   !   ! Leave spatial loops (if any).
   !   _LOOP_END_
   !
   !end subroutine do   

   ! Add model subroutines here.

end module
