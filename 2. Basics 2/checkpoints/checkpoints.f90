! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 2, Exercise 1: Checkpoints
!  Write a subroutine to be used as a checkpoint for debugging.
!  (a) Include a control on a logical variable (Debug=.true. or .false.)
!  (b) Include an additional (optional) string to be printed
!  (c) Include additional (optional) variables to be printed

!> Simple subroutine than prints the file where it is called and the line
!> INPUT: realarg real variable to eventually print
!> OTHER: DEBUG_ if true actually execute the content of the function
!>        line prints the line
 subroutine check_real(DEBUG, realarg, line)
  logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
  real           :: realarg ! optional generic real argument
  integer        :: line

  if (DEBUG .eqv. .TRUE.) then
    print *, 'LINE:',line   ! print file and line
    print*,'arg:', realarg
  end if
end subroutine check_real

#define check_real_(realarg) check_real(DEBUG, realarg,__LINE__)

program test
  logical :: DEBUG = .TRUE.
  real    :: x = 1.23456, y = 6.54321

  ! DEBUG_ is true, this value should be printed
  call check_real_(x)

  DEBUG = .FALSE.
  ! DEBUG_ is now false, this value should NOT be printed
  call check_real_(y)

end program test
