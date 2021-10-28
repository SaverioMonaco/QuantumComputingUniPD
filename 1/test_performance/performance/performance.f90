! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 1, Exercise 3: Test performance
!   Matrix matrix multiplication is many times the bottleneck of linear algebra computations.
!   (a) Write explicitely the matrix-matrix multiplication loop in two different orders.
!   (b) Use the Fortran intrinsic function.
!   (c) Increase the matrix size and use the Fortran Function CPUTIME to monitor the code performance.
!   (d) Use the compiler different optimization flags and monitor the performances

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

end program matrix_mult
