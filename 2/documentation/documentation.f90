! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 2, Exercise 2: Documentation
!  Rewrite Exercise 3 of Week 1 including:
!  (a) Documentation.
!  (b) Comments.
!  (c) Pre- and post- conditions.
!  (d) Error handling
!  (e) Checkpoints

! To compile it:
!   gfortran documentation.f90 -o doc -cpp

!> Matrix Operations
module matrix_utilities
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + matrix_random_initialization:                                      !
  !     Initializes a real matrix randomly                               !
  !     INPUT:                                                           !
  !       > N = integer, number of cols of the matrix                    !
  !       > M = integer, number of rows of the matrix                    !
  !       > rand_range = integer, range of random numbers 0<=x<=range    !
  !     OUTPUT:                                                          !
  !       > mat_A = real*4, dimension(N,M), real random matrix           !
  !                                                                      !
  ! + matrix_multiplication:                                             !
  !     Performs matrix multiplication through loop method               !
  !     INPUT:                                                           !
  !       > mat_A = real*4, dimension(:,:), first input matrix           !
  !       > mat_B = real*4, dimension(:,:), second input matrix          !
  !     OUTPUT:                                                          !
  !       > mat_C = real*4, dimension(:,:), matrix from mat_A * mat_B    !
  !                                                                      !
  ! + matrix_multiplication2:                                            !
  !     Performs matrix multiplication through loop method with the      !
  !     second loop order                                                !
  !     INPUT:                                                           !
  !       > mat_A = real*4, dimension(:,:), first input matrix           !
  !       > mat_B = real*4, dimension(:,:), second input matrix          !
  !     OUTPUT:                                                          !
  !       > mat_C = real*4, dimension(:,:), matrix from mat_A * mat_B    !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES                                                          !
  ! + graphics_printmatrix                                               !
  !   Prints matrix on terminal                                          !
  !   INPUT:                                                             !
  !     > mat_A = real*4, dimension(:,:), matrix to print                !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  implicit none
  contains

  !> Initiale a matrix with random real values
  !> input:  (N, M) size of the matrix
  !> output: A NM-matrix
  function matrix_random_initialization(N, M, rand_range) result(mat_A)
    ! Initialize the matrix
    integer :: N, M, rand_range
    real*4, dimension(N,M) :: mat_A

    ! found this method, probably faster (optimized) than a for loop in a for loop
    call random_number(mat_A)
    mat_A = rand_range*mat_A

  end function matrix_random_initialization

  !> input mat_A, mat_B matrices (2D arrays)
  !> output mat_C array from AB=C. if size not correct returns an empty array
  function matrix_multiplication(mat_A, mat_B) result(mat_C)
    real*4, dimension(:,:)                         :: mat_A, mat_B
    real*4, dimension(size(mat_A,1),size(mat_B,2)) :: mat_C
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
                if (kk == 1) then
                  mat_C(ii,jj) = 0
                end if
                mat_C(ii, jj) = mat_C(ii, jj) + mat_A(ii, kk)*mat_B(kk, jj)
            end do
        end do
      end do
    end if
  end function matrix_multiplication

  !> Other order
  !> input mat_A, mat_B matrices (2D arrays)
  !> output mat_C array from AB=C. if size not correct returns an empty array
  function matrix_multiplication2(mat_A, mat_B) result(mat_C)
    real*4, dimension(:,:)                         :: mat_A, mat_B
    real*4, dimension(size(mat_A,1),size(mat_B,2)) :: mat_C
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
      do jj = 1, size(mat_B,2), 1
        do ii = 1, size(mat_A,1), 1
            do kk = 1, size(mat_B,1), 1
                if (kk == 1) then
                  mat_C(ii,jj) = 0
                end if
                mat_C(ii, jj) = mat_C(ii, jj) + mat_A(ii, kk)*mat_B(kk, jj)
            end do
        end do
      end do
    end if
  end function matrix_multiplication2

  !> Simple function to print a matrix
  subroutine graphics_printmatrix(mat_A)
    integer                 :: ii
    real*4, dimension(:,:)  :: mat_A

    do ii = 1, ubound(mat_A, 1)
      print*, "|", mat_A(ii, :), "|"
    end do

  end subroutine graphics_printmatrix

end module matrix_utilities

!> Module for debugging function
!> Initialize DEBUG_ and set it .TRUE. to be in debug mode
module debugging
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  !  must define a DEBUG_variable that is automatically passed in all    !
  !  debug functions                                                     !
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS: none                                                      !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES                                                          !
  ! + check_real                                                         !
  !   Simple subroutine than prints a real and where it is called and    !
  !   the line                                                           !
  !   INPUT:                                                             !
  !     > realarg = real, real to print                                  !
  !                                                                      !
  ! + check_matrix                                                       !
  !   Simple subroutine than prints matrix and line                      !
  !   INPUT:                                                             !
  !     > mat = real*4, dimension(:,:), matrix to print                  !
  !                                                                      !
  ! + check_custom_matrix_multiplication                                 !
  !   Checks for custom implemented matrix multiplication using matmul   !
  !   INPUT:                                                             !
  !     > mat_A = real*4, dimension(:,:), input matrix for multiplication!
  !     > mat_B = real*4, dimension(:,:), input matrix for multiplication!
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  use matrix_utilities

  contains
  !> Simple subroutine than prints a real and where it is called and the line
  !> INPUT: realarg real variable to eventually print
  !> OTHER: DEBUG_ if true actually execute the content of the function
  !>        line prints the line
  subroutine check_real_(DEBUG, realarg, line)
    logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
    real           :: realarg ! optional generic real argument
    integer        :: line

    if (DEBUG .eqv. .TRUE.) then
      print *, 'LINE:',line   ! print file and line
      print*,'arg:', realarg
    end if
  end subroutine check_real_

  !> Prints matrix and line
  !> INPUT: mat matrix to be printed
  !> OTHER: DEBUG_ if true actually execute the content of the function
  !>        line prints the line
  subroutine check_matrix_(DEBUG, mat, line)
    logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
    integer        :: line
    real*4, dimension(:,:), allocatable :: mat

    if (DEBUG .eqv. .TRUE.) then
      print *, 'LINE:',line   ! print file and line
      call graphics_printmatrix(mat)
    end if
  end subroutine check_matrix_

  !> Computes A times B with my function and compares it to the FORTRAN function
  !> through a subtraction
  !> INPUT: mat_A, mat_B input matrices
  !> OTHER: DEBUG_ if true actually execute the content of the function
  subroutine check_custom_matrix_multiplication_(DEBUG, mat_A, mat_B)
    logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
    real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C, mat_D, mat_E

    ! Initialize matrix according to our function
    mat_C = matrix_multiplication(mat_A,mat_B)

    ! Initialize matrix using fortran function
    mat_D = matmul(mat_A,mat_B)

    mat_E = mat_C - mat_D
    call graphics_printmatrix(mat_E)

  end subroutine check_custom_matrix_multiplication_

end module debugging

! Macro for passing __LINE__ and DEBUG_ automatically
#define check_real(realarg) check_real_(DEBUG, realarg,__LINE__)
#define check_matrix(mat) check_matrix_(DEBUG, mat,__LINE__)
#define check_custom_matrix_multiplication(mat_A, mat_B) check_custom_matrix_multiplication_(DEBUG, mat_A, mat_B)

program matrix_mult
  use debugging
  use matrix_utilities

  real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C, mat_D
  real*8  :: start, finish ! for the CPU times
  logical :: DEBUG = .TRUE.
  integer :: n = 3, m = 4

  allocate(mat_A(n, m))
  allocate(mat_B(m, n))

  mat_A = matrix_random_initialization(n,m,10)
  mat_B = matrix_random_initialization(m,n,10)

  ! TEST call graphics_printmatrix(mat_A)
  call cpu_time(start)  ! start time
  mat_C = matrix_multiplication(mat_A,mat_B)
  call cpu_time(finish) ! end time

  print*
  print*, "Matrix A:"
  call check_matrix(mat_A)

  print*
  print*, "Matrix B:"
  call check_matrix(mat_B)

  print*
  print*, "Matrix C:"
  call check_matrix(mat_C)

  ! Since we are dealing with real values, the matrix will not be exactly0
  ! everywhere but will have small values (even if working properly)
  ! We expect that
  print*
  print*, "Matrix from custom multiplication - Matrix from Fortran multiplication:"
  call check_custom_matrix_multiplication(mat_A, mat_B)

  ! finish - start is Delta T
  print*,'t=', finish - start , "| size=", n      ! print on terminal

  deallocate(mat_A)
  deallocate(mat_B)

  ! TEST for two different orders matrix multiplication
  mat_C = mat_C*0
  print*, ""
  print*, "Test for matrix multiplication in the two different orders"
  mat_A = matrix_random_initialization(n,m,10)
  mat_B = matrix_random_initialization(m,n,10)

  mat_C = matrix_multiplication(mat_A,mat_B)
  mat_D = matrix_multiplication2(mat_A,mat_B)

  print*, "First order A*B"
  call graphics_printmatrix(mat_C)
  print*, "Second order A*B"
  call graphics_printmatrix(mat_D)
  print*, "MATMUL A*B"
  call graphics_printmatrix(matmul(mat_A,mat_B))

end program matrix_mult
