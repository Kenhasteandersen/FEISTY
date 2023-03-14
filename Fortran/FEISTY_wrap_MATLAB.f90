module FEISTY_wrap
  use iso_c_binding
  use FEISTY
  use globals

  implicit none

contains

subroutine f_setupVertical(pprod,nStages,region, bottom, photic) bind(c)

   real(c_double), intent(in), value :: pprod
   real(c_double), intent(in), value :: bottom
   real(c_double), intent(in), value :: photic
   integer(c_int), intent(in), value :: nStages
   integer(c_int), intent(in), value :: region

   call setupVertical(pprod,nStages,region, bottom, photic)

end subroutine f_setupVertical

  subroutine f_calcderivatives(u, dudt) bind(c)
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: dudt(nGrid)

    call calcderivatives(u, dudt)

  end subroutine f_calcderivatives

  end module FEISTY_wrap
