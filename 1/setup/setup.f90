! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 1, Exercise 1: Setup

!   (b) Open emacs and write your first program in FORTRAN
!   (c) Submit a test job.

! Computing pi with a basic montecarlo method
program pi_montecarlo

  implicit none ! types of implicit variable are assumed by their names,
                ! we want to avoid that

  integer :: hits, i          ! hits start from 0, increases by 1 if in circle
  real :: rand_xy(100000,2)   ! random 100000 pair of coordinates

  call random_number(rand_xy) ! All the values of the array are 0<=v<=1 now

  hits = 0
  do i=0,100000,1
   if( rand_xy(i,1)**2 + rand_xy(i,2)**2  < 1) then ! if inside the circle
    hits = hits + 1                                 ! hits increases
   endif
  enddo

  print *, "Pi is approximately:", 4*hits/real(100000)

end program pi_montecarlo

!   (d) (Optional) Connect to the cluster spiro.fisica.unpd.it via ssh
!   and repeat the execution

!   To send the this file through ssh:
!   scp -o KexAlgorithms=diffie-hellman-group1-sha1 ./setup.f90
!   monacos@147.162.55.11:path
