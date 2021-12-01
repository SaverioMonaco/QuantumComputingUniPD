! Having tested matrix multiplication on matrix_mult.f90 I copied the function
! to perform a time analysis depending on the input size

! Operations
module matrix_utilities
  implicit none
  contains
  ! input:  (N, M) size of the matrix
  ! output: A NM-matrix
  function matrix_random_initialization(N, M, rand_range) result(mat_A)
    ! Initialize the matrix
    integer :: N, M, rand_range
    real*4, dimension(N,M) :: mat_A

!    do aa = 1,N,1     ! for each column
!      do bb = 1,M,1   ! for each row
!        A(aa, bb) = RAND(0)*rand_range
!      end do
!    end do

    ! found this method, probably faster (optimized) than a for loop in a for loop
    call random_number(mat_A)
    mat_A = rand_range*mat_A

  end function matrix_random_initialization

  ! input mat_A, mat_B matrices (2D arrays)
  ! output mat_C array from AB=C. if size not correct returns an empty array
  function matrix_multiplication(mat_A, mat_B) result(mat_C)
    real*4, dimension(:,:)                         :: mat_A, mat_B
    real*4, dimension(size(mat_A,2),size(mat_B,1)) :: mat_C
    logical                                        :: check
    integer                                        :: ii,jj,kk

    ! Check if multiplication is possible (shapes)
    if (size(mat_A,2) .eq. size(mat_B,1)) then
      check = .TRUE.
    else
      print*, "Input matrices cannot be multiplied"
      check = .FALSE.
    end if

    ! Initialize output matrix to zeros, when allocating it's not granted every
    ! element is 0
    ! Begin multiplication
    if (check .eqv. .TRUE.) then
      do ii = 1, size(mat_A,1), 1
        do jj = 1, size(mat_B,2), 1
            do kk = 1, size(mat_B,1), 1
                if (kk == 0) then
                  mat_C(ii,jj) = 0
                end if
                mat_C(ii, jj) = mat_C(ii, jj) + mat_A(ii, kk)*mat_B(kk, jj)
            end do
        end do
      end do
    end if
  end function matrix_multiplication

end module matrix_utilities

! Printing and others
module matrix_graphics
  implicit none
  contains

  subroutine graphics_printmatrix(mat_A)
    integer                 :: ii
    real*4, dimension(:,:)  :: mat_A

    do ii = 1, ubound(mat_A, 1)
      print*, "|", mat_A(ii, :), "|"
    end do

  end subroutine graphics_printmatrix

end module matrix_graphics

program matrix_mult

  use matrix_utilities
  use matrix_graphics

  real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C
  real*8 :: start, finish ! for the CPU times

  ! Perform a NxN NxN matrix multiplication
  ! Generate two N by N random matrices
  open(1, file = './loop_flag.csv', status = 'old')
  print*, "t | Size"
  ! We save the times to perform a n by n matrix multiplication
  do n = 5, 100, 1
    allocate(mat_A(n, n))
    allocate(mat_B(n, n))

    mat_A = matrix_random_initialization(n,n,10)
    mat_B = matrix_random_initialization(n,n,10)
    ! TEST call graphics_printmatrix(mat_A)
    call cpu_time(start)  ! start time
    mat_C = matrix_multiplication(mat_A,mat_B)
    call cpu_time(finish) ! end time

    ! finish - start is Delta T
    print*, finish - start , "|", n      ! print on terminal
    write(1,*) finish - start , ",", n   ! save on file

    deallocate(mat_A)
    deallocate(mat_B)

  end do

end program matrix_mult
