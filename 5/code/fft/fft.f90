
module debugmod
! ---------------------------------------------------------------------!
! --------------------------- DOCUMENTATION ---------------------------!
! ---------------------------------------------------------------------!
! TYPES: none                                                          !
! ---------------------------------------------------------------------!
! FUNCTIONS: none                                                      !
! ---------------------------------------------------------------------!
! SUBROUTINES:                                                         !
! + debugging                                                          !
!   INPUT:                                                             !
!   > condition = logic, if condition is TRUE then run the content of  !
!                        this function                                 !
!   > msg = character, optional, message to print                      !
!   > content = some variable to inspect                               !
! ---------------------------------------------------------------------!
! ---------------------------------------------------------------------!
  contains

  subroutine debugging(condition, msg, content)
    logical, intent(IN)                 :: condition
    character(*), intent(IN), optional  :: msg
    class(*), intent(IN), optional      :: content

    if(condition) then
      if (present(content)) then
        select type(content)
          type is (integer(1))
            print*, msg, " => [OK], Variable = ", content
          type is (integer(2))
            print*, msg, " => [OK], Variable = ", content
          type is (integer(4))
            print*, msg, " => [OK], Variable = ", content
          type is (integer(8))
            print*, msg, " => [OK], Variable = ", content
          type is (real(4))
            print*, msg, " => [OK], Variable = ", content
          type is (real(8))
            print*, msg, " => [OK], Variable = ", content
          type is (logical)
            print*, msg, " => [OK], Variable = ", content
        end select
      else
        print*, msg, " => [OK]"
      end if
    end if
    end subroutine

end module debugmod

program fft_test
    use, intrinsic :: iso_c_binding
    use debugmod

    include 'fftw3.f'

    character(20)                   :: folder
    logical                         :: DEBUG
    integer                         :: signal_size,iostat,ii

    integer*8                         :: plan
    double complex,   dimension(:), allocatable   :: signal_out, signal
    double precision, dimension(:), allocatable   :: real_signal

    DEBUG = .TRUE.

    print*, "--------------------------------------------"
    print*, "+          FOURIER TRANSFORM TEST          +"
    print*, "--------------------------------------------"

    write(*,"(A)",advance='no') " + Select file: "
    read (*,*) folder

    if(DEBUG) then
      print*, "+ Folder:", folder
    endif

    ! DATA LOADING
    ! Reading the file line by line to get the size of the instances
    open(10, file = folder)
    signal_size = 0
    do while(1.eq.1)
         read(10,*,iostat = iostat)
         if (iostat /= 0) then
           exit
         end if
         signal_size = signal_size + 1
    end do
    close(10)

    print*, "--------------------------------------------"
    print*, "+ Loading Data"
    print*, "+   Size of Data:", signal_size

    allocate(real_signal(signal_size))
    allocate(signal(signal_size))

    open(11, file = folder)
    ! then we use the signal_size found and use it to fill the signal_array
    do ii=1, signal_size, 1
      read(11,*) real_signal(ii)
    end do
    close(11)
    signal = cmplx(real_signal,0*real_signal,kind=8)

    print*, "+ Data loaded"
    print*, "+"

    if(DEBUG) then
      print*, "+ Signal:"
      do ii=1, 5, 1
        print*, "+", signal(ii)
      end do
      print*, "+     ..."
    end if

    print*, "+"
    print*, "+ Performing Fourier Transform..."

    allocate(signal_out(signal_size))

    call dfftw_plan_dft_1d(plan, signal_size, signal, signal_out, FFTW_FORWARD, FFTW_ESTIMATE)
    call debugging(DEBUG, '+   1/2')

    call dfftw_execute(plan)!, signal, signal_out)
    call debugging(DEBUG, '+   2/2')

    print*, "+ DONE!"
    print*, "+"

    if(DEBUG) then
      print*, "+ Output Signal:"
      do ii=1, 5, 1
        print*, "+", signal_out(ii)
      end do
      print*, "+     ..."
    end if

    print*, "+"
    write(*,"(A)",advance='no') " + Writing on file"

    open(14, file = 'output_signal.csv')
    do ii = 1, size(signal_out), 1
      write(14, *) signal_out(ii)
    end do
    close(14)
    print*, "DONE!"

    call dfftw_destroy_plan(plan)

    print*, "+"
    print*, "+-------------------------------------------"
    print*, "+-------------------------------------------"
end program fft_test
