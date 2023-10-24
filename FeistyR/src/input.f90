!
! Module to handle the reading of the input file
!
module input
  implicit none
  character(256) :: file_path, file_path_V


  contains


  subroutine open_inputfile(file_unit, io_err)
    integer,  intent(out) :: file_unit, io_err
    open(newunit=file_unit,action='read', file=file_path,iostat=io_err)!C:/Users/Admin/Desktop/FEISTY-main/FEISTY-main../FeistyR/data
    call check_iostat(io_err, &
                      "Could not open file 'input.nml', perhaps it does not exist?")
  end subroutine open_inputfile

  subroutine close_inputfile(file_unit, io_err)
    integer,  intent(in) :: file_unit, io_err
    call check_iostat(io_err, &
                      "Error reading name-list, please correct file content")
    close (file_unit)
  end subroutine close_inputfile

    ! -----------------------------------------------
    ! Find problems with input
    ! -----------------------------------------------

  subroutine check_iostat(iostat, msg)
    integer, intent(in) :: iostat
    character(len=*), intent(in) :: msg

    if ( iostat == 0 ) return

    write(*,'(a,i0,/,tr3,a)') 'ERROR = ', iostat, msg
    stop

  end subroutine check_iostat



end module input

! 
! Routines used to pass the paths for the input files from R to the Fortran library:
!
    subroutine passpath(length, file_path_in) bind(C)
      use input
     use iso_c_binding, only: c_int, c_char
     integer(c_int), intent(in):: length ! Length of the character string being sent in
     character(c_char), dimension(length), intent(in) :: file_path_in
     integer:: i

    ! Transfor from the R-file_path_in to the Fortran format
     do i = 1, length
       file_path(i:i) = file_path_in(i)
     end do
    end subroutine passpath

    subroutine passpathv(length, file_path_in) bind(C)
      use input
      use iso_c_binding, only: c_int, c_char
      integer(c_int), intent(in):: length ! Length of the character string being sent in
      character(c_char), dimension(length), intent(in) :: file_path_in
      integer:: i
 
     ! Transfor from the R-file_path_in to the Fortran format
      do i = 1, length
        file_path_V(i:i) = file_path_in(i)
      end do
    end subroutine passpathv

