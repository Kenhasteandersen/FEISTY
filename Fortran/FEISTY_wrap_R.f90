! for R
subroutine f_setupbasic(szprod,lzprod,bprod, Ts, Tb)
   use FEISTY !, only:
   use globals

   real(dp), intent(in)::szprod,lzprod,bprod, Ts, Tb

   call setupbasic(szprod,lzprod,bprod, Ts, Tb)
end subroutine f_setupbasic

subroutine f_setupbasic2(szprod,lzprod, bprod, nStages, Ts, Tb,etaMature)
   use FEISTY ! only : setupbasic2
   use globals

   real(dp), intent(in):: szprod,lzprod, bprod, Ts, Tb, etaMature
   integer, intent(in):: nStages

   call setupbasic2(szprod,lzprod, bprod, nStages, Ts, Tb, etaMature)
end subroutine f_setupbasic2

subroutine f_setupvertical(szprod,lzprod, bent, nStages, region, bottom, photic, etaMature)
   use FEISTY ! only :
   use globals

   real(dp), intent(in):: szprod,lzprod, bottom, photic, bent, etaMature
   integer, intent(in):: nStages, region

   call setupvertical(szprod,lzprod, bent, nStages, region, bottom, photic, etaMature)
end subroutine f_setupvertical

subroutine f_setupverticalglobal(szprod, lzprod, bprod, bottom, photic, dgrid, tprof, nStages)
   use FEISTY ! only :
   use globals

   real(dp), intent(in):: szprod, lzprod, bprod, bottom, photic, dgrid(200), tprof(200)
   integer, intent(in):: nStages

   call setupVerticalGlobal(szprod, lzprod, bprod, bottom, photic, Dgrid, Tprof, nStages)
end subroutine f_setupverticalglobal

subroutine f_setupsquid(szprod,lzprod, bottom, nStages)
   use FEISTY ! only :
   use globals

   real(dp), intent(in):: szprod,lzprod, bottom
   integer, intent(in) :: nStages

   call setupsquid(szprod,lzprod, bottom, nStages)
end subroutine f_setupsquid

subroutine f_calcderivatives(u, dudt)
   use FEISTY !, only:
   use globals

   real(dp), intent(in):: u(nGrid)
   real(dp), intent(out):: dudt(nGrid)

   call calcderivatives(u, dudt)
end subroutine f_calcderivatives

subroutine f_calcderivativesSquid(u, dudt)
   use FEISTY !, only:
   use globals

   real(dp), intent(in):: u(nGrid)
   real(dp), intent(out):: dudt(nGrid)

   call calcderivativesSquid(u, dudt)
end subroutine f_calcderivativesSquid

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

subroutine f_simulateEuler(u, dudt, tEnd, dt)
    use FEISTY
    use globals

         real(dp),intent(inout) :: u(nGrid)
         real(dp),intent(in):: dudt(nGrid)
         real(dp),intent(in) :: tEnd, dt

         call simulateEuler(u, dudt, tEnd, dt)

end subroutine f_simulateEuler
