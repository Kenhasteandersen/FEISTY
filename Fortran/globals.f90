!
! from NUM library
! subroutines of temperature may be used for further work
! new added useful functions : linspace, diag
!

module globals
   use input
   implicit none
   integer, parameter :: dp = kind(0.d0) ! double precision
   !
   ! Useful mathematical constants:
   !
   real(dp), parameter :: onethird = 1.d0/3.d0
   real(dp), parameter :: twothirds = 2.d0/3.d0
   real(dp), parameter :: threequarters = 3.d0/4.d0
   real(dp), parameter :: pi = 4*ATAN(1.d0)

   ! Small number to avoid divisions by zero
   real(dp), parameter :: eps = 1d-200

   ! Temperature Q10 corrections
   real(dp), parameter :: Q10=1.88d0, Q10m=1.88d0, Q10mPetrik=2.35d0
   real(dp) :: fTemp, fTempm, fTempdem, fTempmdem
   real(dp), parameter:: Tref = 10.d0 ! Reference temperature

   real(dp), allocatable :: fTempV(:), fTempmV(:)

   !
   ! Specification of what to do with HTL losses:
   !
   real(dp) :: fracHTL_to_N ! Half becomes urine that is routed back to N
   real(dp) :: fracHTL_to_POM ! Another half is fecal pellets that are routed back to the largest POM size class
   real(dp) :: rhoCN
   !real(dp), parameter :: fracHTL_to_N = 0.5 ! Half becomes urine that is routed back to N
   !real(dp), parameter :: fracHTL_to_POM = 0.5 ! Another half is fecal pellets that are routed back to the largest POM size class
   !real(dp), parameter :: rhoCN = 5.68

contains

   ! -----------------------------------------------
   ! Read in general parameters
   ! -----------------------------------------------
   subroutine read_namelist_general()
      integer :: file_unit, io_err

      namelist /input_general/ rhoCN, fracHTL_to_N, fracHTL_to_POM

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_general, iostat=io_err)
      call close_inputfile(file_unit, io_err)

   end subroutine read_namelist_general

   ! -----------------------------------------------
   ! Temperature Q10 function
   ! -----------------------------------------------
   function calfTemp(Q10, T) result(f) !calculate temperature factor
      real(dp), intent(in):: Q10, T
      real(dp):: f

      f = Q10**((T - Tref)/10.d0)
   end function calfTemp

   ! -----------------------------------------------
   ! Update the temperature corrections only if T has changed
   ! -----------------------------------------------
!   subroutine updateTemperature(T)
!      real(dp), intent(in) :: T
!      real(dp), save :: Told = -1000.
!
!      if (T .ne. Told) then
!         Told = T
!         fTemp2 = fTemp(2.d0, T)
!         fTemp15 = fTemp(1.5d0, T)
!      end if
!   end subroutine updateTemperature

! update temperature setupbasic (Petrik et al., 2019) and setupbasic2
subroutine updateTemp(Ts, Tb)
      real(dp), intent(in) :: Ts, Tb
      real(dp), save :: Tolds = -1000.d0
      real(dp), save :: Toldb = -1000.d0

      if (Ts .ne. Tolds .OR. Tb .ne. Toldb) then
         Tolds = Ts
         Toldb = Tb
         fTemp = calfTemp(Q10, Ts)  !Q10=1.88 clearance rate,  maximum consumption rate
         fTempm = calfTemp(Q10mPetrik, Ts) !Q10m=2.35 for metabolism    Petrik
         fTempdem = calfTemp(Q10, Tb)  !for demersal
         fTempmdem = calfTemp(Q10mPetrik, Tb) !for demersal

      end if

end subroutine updateTemp

! update Temperature for vertical version van Denderen et al., 2021
subroutine updateTempV(depthDay, depthNight, bottom, region)
 real(dp), intent(in) :: depthDay(:, :), depthNight(:, :), bottom
 integer :: i,region
 real(dp), allocatable :: dist(:,:), TQ10(:), TQ10m(:), fTemp_stepV(:,:), fTempm_stepV(:,:)
 real(dp) :: tempdata(5501,5) ! contain tempdata from van Denderen et al., 2021 + default(10 celcius)

    if (allocated (fTempV)) then
        deallocate (fTempV)
        deallocate (fTempmV)
    end if

allocate (dist(size(depthDay,1), size(depthDay,2)))
allocate (TQ10(int(bottom+1)))
allocate (TQ10m(int(bottom+1)))
allocate (fTemp_stepV(size(depthDay,1), size(depthDay,2)))
allocate (fTempm_stepV(size(depthDay,1), size(depthDay,2)))
allocate (fTempV(size(depthDay,2)))
allocate (fTempmV(size(depthDay,2)))

    dist = 0.d0
    TQ10 = 0.d0
    TQ10m = 0.d0
    fTemp_stepV = 0.d0
    fTempm_stepV = 0.d0
    fTempV = 0.d0
    fTempmV  = 0.d0

open(unit=1,action='read', file="../input/tempdata.dat",status="old")!C:/Users/Admin/Desktop/FEISTY-main/FEISTY-main
do i = 1,5501
    read(1,*) tempdata(i,1),tempdata(i,2),tempdata(i,3),tempdata(i,4) !depth 0-5500(no use), tropical, temperate, boreal, default(10 celcius)
end do
close(1)
tempdata(:,5)=10.d0 !default temp, so no temp-effects

dist = (depthDay + depthNight)/2.d0
! region+1: 1+1 tropical, 2+1 temperate, 3+1 boreal, 4+1 default(10 celcius)
TQ10 =  Q10**((tempdata(1:bottom+1 , (region+1))-10.d0)/10.d0)
TQ10m =  Q10m**((tempdata(1:bottom+1 , (region+1))-10.d0)/10.d0)


do i=1,size(dist,2)
fTemp_stepV(:,i) = dist(:,i) * TQ10
fTempm_stepV(:,i) = dist(:,i) * TQ10m
end do

fTempV = sum(fTemp_stepV,1)
fTempmV = sum(fTempm_stepV,1)

end subroutine updateTempV

! update temperature for global offline coupling (vertical version)
! Dgrid: Depth grids,from earth system models or dataset, from surface to bottom (boundary data, not grid center data)
! Tprof: Temperature profile, from earth system models or dataset, from surface T to bottom T (boundary data, not grid center data)
! bottom: Seafloor depth (boundary)
subroutine updateTempVG(Dgrid, Tprof, depthDay, depthNight, bottom)
      real(dp), intent(in) :: Dgrid(:), Tprof(:),depthDay(:, :),depthNight(:, :),bottom
      real(dp), allocatable :: tempdata (:),temp(:), dist(:,:), TQ10(:), TQ10m(:), fTemp_stepV(:,:), fTempm_stepV(:,:)
      integer :: i, idxn
      real(dp), save :: Told = -1000.d0

! if temp changes (location changes) then update
! if (sum(Tprof) .ne. Told) then
!         Told = sum(Tprof)

    if (allocated (tempdata)) then
        deallocate (tempdata)
        deallocate (dist)
        deallocate (TQ10)
        deallocate (TQ10m)
        deallocate (fTemp_stepV)
        deallocate (fTempm_stepV)
        deallocate (fTempV)
        deallocate (fTempmV)
    end if

    if (allocated (fTempV)) then
        deallocate (fTempV)
        deallocate (fTempmV)
    end if

    allocate (tempdata(int(bottom+1)))
    allocate (dist(size(depthDay,1), size(depthDay,2)))
    allocate (TQ10(int(bottom+1)))
    allocate (TQ10m(int(bottom+1)))
    allocate (fTemp_stepV(size(depthDay,1), size(depthDay,2)))
    allocate (fTempm_stepV(size(depthDay,1), size(depthDay,2)))
    allocate (fTempV(size(depthDay,2)))
    allocate (fTempmV(size(depthDay,2)))

    tempdata = 0.d0
    dist = 0.d0
    TQ10 = 0.d0
    TQ10m = 0.d0
    fTemp_stepV = 0.d0
    fTempm_stepV = 0.d0
    fTempV = 0.d0
    fTempmV  = 0.d0

idxn = minloc(abs(Dgrid - bottom),dim=1) ! find index that corresponds to seafloor

do i = 1,idxn-1
    if (allocated(temp)) then
          deallocate(temp)
    end if
  allocate (temp(int(Dgrid(i+1)-Dgrid(i)+1)))
    temp = linspace(Tprof(i),Tprof(i+1),int((Dgrid(i+1)-Dgrid(i)+1)))
    tempdata(int(Dgrid(i)+1):int(Dgrid(i+1))) = temp(1:(size(temp)-1))

end do
tempdata(size(tempdata)) = Tprof(idxn)

dist = (depthDay + depthNight)/2.d0
TQ10 =  Q10**((tempdata(1:int(bottom+1))-10.d0)/10.d0)
TQ10m =  Q10m**((tempdata(1:int(bottom+1))-10.d0)/10.d0)

do i=1,size(dist,2)
fTemp_stepV(:,i) = dist(:,i) * TQ10
fTempm_stepV(:,i) = dist(:,i) * TQ10m
end do

fTempV = sum(fTemp_stepV,1)
fTempmV = sum(fTempm_stepV,1)

!end if

end subroutine updateTempVG


  ! linspace function:
      function linspace(value_start, value_end, length) result(res)
         real(dp), intent(in)::value_start, value_end
         integer, intent(in)::length
         real(dp) ::dx, res(length)
         integer::i

         dx = (value_end - value_start)/(length - 1)
         res(1:length) = [(value_start + (i - 1)*dx, i=1, length)]
      end function linspace

      ! a simple function to create diagonal matrix, only for 1-D array
    function diag(A) result(res)
    real(dp), intent(in) :: A(:)
    real(dp) :: res(size(A),size(A))
    integer :: i
    res=0.d0
    do i = 1, size(A)
        res(i,i) = A(i)
    end do
    end function

end module globals
