module FEISTY_wrap
  use iso_c_binding
  use FEISTY
  use globals

  implicit none

contains

subroutine f_setupVertical(pprod,nStages,region) bind(c)

   real(c_double), intent(in), value :: pprod
   integer(c_int), intent(in), value :: nStages
   integer(c_int), intent(in), value :: region

   call setupVertical(pprod,nStages,region)

end subroutine f_setupVertical

  subroutine f_calcderivatives(u, dudt) bind(c)
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: dudt(nGrid)

    call calcderivatives(u, dudt)

  end subroutine f_calcderivatives

  end module FEISTY_wrap
