! for R
subroutine f_setupFEISTY(pprod)
   use FEISTY !, only:setupFEISTY
   use globals

   real(dp), intent(in)::pprod

   call setupFEISTY(pprod)
end subroutine f_setupFEISTY

subroutine f_calcderivatives(u, dudt)
   use FEISTY !, only: calcderivatives, nGrid
   use globals

   real(dp), intent(in):: u(nGrid)
   real(dp), intent(out):: dudt(nGrid)

   call calcderivatives(u, dudt)
end subroutine f_calcderivatives

subroutine f_getrates(u, dudt, flvl_r, mortpred_r, g_r)
use FEISTY
use globals

         real(dp),intent(in) :: u(nGrid)
         real(dp),intent(inout):: dudt(nGrid)
         real(dp),intent(out) :: flvl_r(nGrid)
         real(dp),intent(out) :: mortpred_r(nGrid)
         real(dp),intent(out) :: g_r(nGrid-nResources)

         call getrates(u, dudt, flvl_r, mortpred_r, g_r)

end subroutine f_getrates

