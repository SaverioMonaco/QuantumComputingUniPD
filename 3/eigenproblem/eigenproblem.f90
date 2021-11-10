! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 3, Exercise 2: Eigenproblem
! Consider a random Hermitian matrix A of size N .
!   (a) Diagonalize A and store the N eigenvalues λ i in crescent order.
!   (b) Compute the normalized spacings between eigenvalues

! Hint: Use LAPACK (ZHEEV)

module cmatrices
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES:                                                               !
  ! + cmatrix:                                                           !
  !              Structure for complex 2D matrices                       !
  !              > dim:     dimension of the matrix                      !
  !              > element: actual array of elements                     !
  !              > trace:   trace of the matrix                          !
  !              > det:     determinant                                  !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + cmatrix_trace:                                                     !
  !     Computes the trace                                               !
  !     INPUT:                                                           !
  !       > cmat = type(cmatrix), input complex matrix                   !
  !     OUTPUT:                                                          !
  !       > trace = complex*16, trace of thecomplex matrix               !
  !                                                                      !
  ! + cmatrix_randinit:                                                  !
  !     Initialiazes complex matrix randomly                             !
  !     INPUT:                                                           !
  !       > nrow = integer, number of rows of matrix                     !
  !       > ncol = integer, number of columns of matrix                  !
  !       > range = real, range of numbers: 0<=x<=range                  !
  !     OUTPUT:                                                          !
  !       > cmat = type(cmatrix), complex matrix                         !
  !                                                                      !
  ! + cmatrix_randinit_hermitian:                                        !
  !     Initialiazes complex hermitian matrix randomly                   !
  !     INPUT:                                                           !
  !       > n = integer, dimension of matrix                             !
  !       > range = real, range of numbers: 0<=x<=range                  !
  !     OUTPUT:                                                          !
  !       > cmat = type(cmatrix), complex hermitian matrix               !
  !                                                                      !
  ! + cmatrix_init:                                                      !
  !     Initialiazes complex matrix type given the 2d array              !
  !     INPUT:                                                           !
  !       > array2d = real*16, dimension(:,:), input 2d array            !
  !     OUTPUT:                                                          !
  !       > cmat = type(cmatrix), complex matrix                         !
  !                                                                      !
  ! + cmatrix_adjoint                                                    !
  !     Outputs the adjoints of the input complex matrix                 !
  !     INPUT:                                                           !
  !       > cmat = type(cmatrix), input complex matrix                   !
  !     OUTPUT:                                                          !
  !       > cmat = type(cmatrix), adjoint matrix                         !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES                                                          !
  ! + cmatrix_print                                                      !
  !   Prints matrix on terminal                                          !
  !   INPUT:                                                             !
  !     > cmat = type(cmatrix), matrix to print                          !
  !                                                                      !
  ! + cmatrix_write                                                      !
  !   Writes matrix to file                                              !
  !   INPUT:                                                             !
  !     > cmat = type(cmatrix), matrix to store                          !
  ! ---------------------------------------------------------------------!
  implicit none

  !> Here we define the type cmatrix, a double complex 2D matrix of generic
  !> dimension
  type cmatrix
    integer, dimension(2)                   :: dim        ! dimension of the matrix
    complex*16, dimension(:,:), allocatable :: element
    complex*16                              :: trace, det
  end type cmatrix

  interface operator(.Adj.)
    module procedure cmatrix_adjoint
  end interface

  interface operator(.Trace.)
    module procedure cmatrix_trace
  end interface

  interface operator(.Eigens.)
    module procedure cmatrix_eigenvalues
  end interface

  interface operator(.Eigenspacing.)
    module procedure cmatrix_eigenspacing
  end interface

  contains

  function cmatrix_trace(cmat) result(trace)
    integer                   :: ii
    complex*16                :: trace
    type(cmatrix), intent(IN) :: cmat

    if(cmat%dim(1) == cmat%dim(2)) then ! iff the matrix is square
      trace = 0
      do ii = 1, size(cmat%element,1), 1
        trace = trace + cmat%element(ii,ii)
      end do
    else
      ! Weird way to assign trace as NaN
      trace = 0
      trace = trace/trace
    end if
  end function cmatrix_trace

  function cmatrix_randinit(nrow,ncol,range) result(cmat)
    integer                 :: nrow, ncol
    real                    :: range
    real*16, dimension(:,:), allocatable :: elem_real,elem_imag

    type(cmatrix)          :: cmat

    cmat%dim = (/ nrow, ncol /) ! Assign the dimension

    ! To define the complex matrix, we initialize 2 real matrices for the
    ! real and imaginary components and then we combine them together
    allocate(elem_real(ncol,nrow))
    allocate(elem_imag(ncol,nrow))

    call random_number(elem_real)
    call random_number(elem_imag)

    elem_real = range*elem_real
    elem_imag = range*elem_imag

    ! Combining the matrices
    cmat%element = cmplx(elem_real,elem_imag)
    cmat%trace   = cmatrix_trace(cmat)
    cmat%det     = 0
  end function cmatrix_randinit

  function cmatrix_randinit_hermitian(n,range) result(cmat)
    integer                 :: n
    real                    :: range
    real*16, dimension(:,:), allocatable :: elem_real,elem_imag
    integer                 :: ii,jj,kk

    type(cmatrix)          :: cmat

    cmat%dim = (/ n, n /) ! Assign the dimension

    ! To define the complex matrix, we initialize 2 real matrices for the
    ! real and imaginary components and then we combine them together
    allocate(elem_real(n,n))
    allocate(elem_imag(n,n))

    ! Imaginary part diagonal
    do ii = 1, n, 1
      call random_number(elem_real(ii, ii))
      call random_number(elem_imag(ii, ii))
    end do

    ! Upper right part
    do ii = 1, n, 1
      do jj = ii + 1, n, 1
        call random_number(elem_real(jj, ii))
        call random_number(elem_imag(jj, ii))
      end do
    end do

    ! Flip negatively imaginary values
    do ii = 2, n, 1
      do jj = 1,ii-1, 1
        elem_real(jj, ii) = + elem_real(ii,jj)
        elem_imag(jj, ii) = - elem_imag(ii,jj)
      end do
    end do

    elem_real = range*elem_real
    elem_imag = range*elem_imag

    ! Combining the matrices
    cmat%element = cmplx(elem_real,elem_imag)
    cmat%trace   = cmatrix_trace(cmat)
    cmat%det     = 0
  end function cmatrix_randinit_hermitian

  function cmatrix_init(array2d) result(cmat)
    integer                 :: nrow, ncol
    real                    :: range
    real*16, dimension(:,:), allocatable :: elem_real,elem_imag
    complex*16, dimension(:,:)           :: array2d

    type(cmatrix)          :: cmat

    cmat%dim = (/ size(array2d,1), size(array2d,2) /) ! Assign the dimension

    cmat%element = array2d
    cmat%trace   = cmatrix_trace(cmat)
    cmat%det   = 0
  end function cmatrix_init

  function cmatrix_adjoint(cmat) result(cmat_adj)
    type(cmatrix), intent(IN) :: cmat
    type(cmatrix)             :: cmat_adj

    cmat_adj = cmatrix_init( conjg(transpose(cmat%element)) )
  end function cmatrix_adjoint

  function cmatrix_eigenvalues(cmat) result(eigens)
    type(cmatrix), intent(IN)                   :: cmat
    double precision, dimension(:), allocatable :: eigens, RWORK
    integer                                     :: INFO, LWORK
    integer                                     :: N
    integer, parameter                          :: LWMAX = 100000
    complex*16                                  :: WORK(LWMAX)

    ! Check if matrix is squared
    if(cmat%dim(1) == cmat%dim(2)) then
      N = cmat%dim(1)

      allocate(eigens(N))
      allocate(RWORK(3*N-2))

      ! Compute optimal size of workspace
      LWORK = -1
      call ZHEEV('N', 'U', N, cmat%element, N, eigens, WORK, LWORK, RWORK, INFO)
      LWORK = min(LWMAX, int(WORK(1)))

      ! Compute eigenvalues
      call ZHEEV('N', 'U', N, cmat%element, N, eigens, WORK, LWORK, RWORK, INFO)
    end if
  end function cmatrix_eigenvalues

  function cmatrix_eigenspacing(cmat) result(spacing)
    type(cmatrix), intent(IN)                   :: cmat
    double precision, dimension(:), allocatable :: eigens, spacing
    double precision                            :: inverseaveragelambda
    integer                                     :: ii, N

    ! Check if matrix is squared
    if(cmat%dim(1) == cmat%dim(2)) then
      N = cmat%dim(1)
      allocate(eigens(N))
      allocate(spacing(N-1))

      eigens = cmatrix_eigenvalues(cmat)
      inverseaveragelambda = (N-1)/ (eigens(N)-eigens(1))
      
      do ii = 1, N-1, 1
        spacing(ii) = inverseaveragelambda * (eigens(ii+1) - eigens(ii))
      end do

    end if
  end function cmatrix_eigenspacing

  subroutine cmatrix_print(cmat)
    type(cmatrix) :: cmat
    integer       :: ii

    do ii = 1, ubound(cmat%element, 1)
      print*, "|", cmat%element(ii, :), "|"
    end do

  end subroutine cmatrix_print

  subroutine cmatrix_write(filename, cmat)
    character(len = 25) :: filename
    type(cmatrix)       :: cmat
    integer             :: ii

    open(1, file = filename, status = 'old')

    do ii = 1, ubound(cmat%element, 1)
      write(1, '(*(F0.4, 1x, SP, F0.4, "i", 2x))')  cmat%element(ii, :)
    end do
  end subroutine cmatrix_write
end module cmatrices

program eigenvalue
  use cmatrices

  type(cmatrix) :: cmat

  cmat  = cmatrix_randinit_hermitian(3,10.0)
  call cmatrix_print(cmat)

  ! (a) Diagonalize A and store the N eigenvalues λ i in crescent order.
  print*, "Eigenvalues (ascending):", .Eigens.cmat

  ! (b) Compute the normalized spacings between eigenvalues
  print*, "Normalized spacings:", .Eigenspacing.cmat


end program eigenvalue
