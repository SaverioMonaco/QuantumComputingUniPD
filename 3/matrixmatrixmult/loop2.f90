! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 3, Exercise 1: Scaling of the matrix-matrix multiplication
! Consider the program developed in Exercise 3 of Week 1 (matrix-matrix multiplication).
!   (a) Write a python script that changes N between two values N min and N max , and launches the
!       program.
!   (b) Store the results of the time needed in different files depending on the multiplication method used.
!   (c) Fit the scaling of the time needed for different methods as a function of the input size. Consider
!       the biggest possible difference between N min and N max .
!   (d) Plot the results for the different multiplication methods.

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

program matrix_mult

  use matrix_utilities

  real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C
  real*8                              :: start, finish ! for the CPU times
  integer                             :: Nstart, Nfinish, Ndelta

  open(1, file = './loop2.csv', status = 'old')
  print *, "Enter N start, N finish and Delta N"
  read (*,*) Nstart,Nfinish,Ndelta

  ! We save the times to perform a n by n matrix multiplication
  do n = Nstart, Nfinish, Ndelta
    allocate(mat_A(n, n))
    allocate(mat_B(n, n))

    mat_A = matrix_random_initialization(n,n,10)
    mat_B = matrix_random_initialization(n,n,10)
    ! TEST call graphics_printmatrix(mat_A)
    call cpu_time(start)  ! start time
    mat_C = matrix_multiplication2(mat_A,mat_B)
    call cpu_time(finish) ! end time

    ! finish - start is Delta T
    print*, finish - start , "|", n      ! print on terminal
    write(1,*) finish - start , ",", n   ! save on file

    deallocate(mat_A)
    deallocate(mat_B)

  end do

end program matrix_mult
