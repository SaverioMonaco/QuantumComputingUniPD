! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 1, Exercise 2: Number precision
!   Integer and real numbers have a finite precision.
!   Explore the limits of INTEGER and REAL in Fortran.
!   Sum the numbers 2.000.000 and 1 with INTEGER*2 and INTEGER*4

program n_prec

  implicit none

  integer*2 :: big_number   ! 2000000
  integer*4 :: small_number ! 1

  integer*4 :: result
  integer*2 :: maximum_value, maximum_value_plus1

  big_number   = 2000000
  small_number = 1

  result = big_number + small_number
  print *, "Adding 2000000 + 1 gives", result

  maximum_value = 32767
  print *, "The maximum number for a INTEGER*2 is", maximum_value

  maximum_value_plus1 = maximum_value + 1
  print *, "When adding just by 1", maximum_value, "becomes", maximum_value_plus1, "<- OVERFLOW"
end program n_prec

! an error appears when compiling.
! in order to compile it anyway you must add the flag -fno-range-check:
! gfortran 2a.f90 -o 2a -fno-range-check
