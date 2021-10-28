! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 1, Exercise 2: Number precision
!   Integer and real numbers have a finite precision.
!   Explore the limits of INTEGER and REAL in Fortran.
!   Sum the numbers pi·10^32 and 2·10^21 in single and double precision.

program n_prec2

  implicit none

  ! PI = 4.D0*DATAN(1.D0)
  ! This style ensures that the maximum precision available on ANY architecture
  ! is used when assigning a value to PI. - StackOverflow

  real                        :: single = 1.12345678901234567890
  double precision            :: double = 1.12345678901234567890

  real                        :: real_big_pi = 4.D0*DATAN(1.D0)*10**32
  real                        :: real_big_sqrt2 = sqrt(real(2))*10**21

  double precision            :: double_big_pi = 4.D0*DATAN(1.D0)*10**32
  double precision            :: double_big_sqrt2 = sqrt(real(2))*10**21

  print *, "Single precision:   ", single
  print *, "Double precision:   ", double

  print *, "Real multiplication:", real_big_pi*real_big_sqrt2
  print *, "Real multiplication:", double_big_pi*double_big_sqrt2

end program n_prec2
