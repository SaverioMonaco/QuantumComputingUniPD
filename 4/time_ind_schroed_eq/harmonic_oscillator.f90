! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 4: Solving time independent Schroedinger Equaton
! To compile it:
!   gfortran eigenproblem.f90 -o eigen -llapack

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
  !                                                                      !
  ! + cmatrix_heigenvalues                                               !
  !     Creates an array of the eigenvalues of a complex hermitian matrix!
  !     INPUT:                                                           !
  !       > cmat = type(cmatrix), input complex matrix  (must be herm)   !
  !     OUTPUT:                                                          !
  !       > eigens = double precision, array of eigenvalues              !
  !                                                                      !
  ! + cmatrix_heigenspaces                                               !
  !     Creates an array of the normalized spacings of eigenvalues       !
   !    of a complex hermitian matrix                                    !
  !     INPUT:                                                           !
  !       > cmat = type(cmatrix), input complex matrix  (must be herm)   !
  !     OUTPUT:                                                          !
  !       > spacing = double precision, array of eigenvalues             !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES                                                          !
  ! + cmatrix_heigens                                                    !
  !   Computes both eigenvalues and eigenvectors                         !
  !   INPUT:                                                             !
  !     > cmat = type(cmatrix), input matrix                             !
  !     > eigenv = complex(kind=8), dimension(:), eigenvalues OVERWR     !
  !     > eigenh = complex(kind=8), dimension(:,:), eigenvectors OVERWR  !
  !                                                                      !
  ! + cmatrix_print                                                      !
  !   Prints matrix on terminal                                          !
  !   INPUT:                                                             !
  !     > cmat = type(cmatrix), matrix to print                          !
  !                                                                      !
  ! + cmatrix_write                                                      !
  !   Writes matrix to file                                              !
  !   INPUT:                                                             !
  !     > cmat = type(cmatrix), matrix to store                          !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  implicit none

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

  subroutine cmatrix_herm_eigens(cmat,eigenv,eigenh,success)
    type(cmatrix)                               :: cmat
    real*8, dimension(:)                        :: eigenv
    complex(kind=8), dimension(:,:)             :: eigenh
    integer, optional                           :: success

    ! LAPACK variables
    double precision, dimension(:), allocatable   :: RWORK
    integer                                       :: INFO, LWORK
    integer                                       :: N
    integer, parameter                            :: LWMAX = 100000
    complex*16                                    :: WORK(LWMAX)
    complex(kind=8), dimension(:,:), allocatable  :: VR
    ! Check if matrix is squared
    if(cmat%dim(1) == cmat%dim(2)) then
      N = cmat%dim(1)

      allocate(RWORK(3*N-2))
      allocate(VR(N,N))

      ! Compute optimal size of workspace
      LWORK = -1
      eigenh = cmat%element

      call ZHEEV('Vectors', 'U', N, eigenh, N, eigenv, WORK,LWORK,RWORK,INFO)
      LWORK = min(LWMAX, int(WORK(1)))

      ! Compute eigenvalues
      call ZHEEV('Vectors', 'U', N, eigenh, N, eigenv, WORK,LWORK,RWORK,INFO)

      if(present(success)) then
        success = INFO
      end if
    end if
  end subroutine cmatrix_herm_eigens

  subroutine cmatrix_print(cmat)
    type(cmatrix) :: cmat
    integer       :: ii

    do ii = 1, ubound(cmat%element, 1)
      print*, "|", cmat%element(ii, :), "|"
    end do

  end subroutine cmatrix_print

  subroutine cmatrix_print_real(cmat)
    type(cmatrix) :: cmat
    integer       :: ii

    do ii = 1, ubound(cmat%element, 1)
      print*, "|", real(cmat%element(ii, :)), "|"
    end do

  end subroutine cmatrix_print_real

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

module qho
! ---------------------------------------------------------------------!
! --------------------------- DOCUMENTATION ---------------------------!
! ---------------------------------------------------------------------!
! TYPES: none                                                          !
! ---------------------------------------------------------------------!
! FUNCTIONS:                                                           !
!  + qho_H_init:                                                       !
!    Initializes the discreet hamiltonian matrix for a quantum         !
!    harmonic oscillator (QHO)                                         !
!    INPUT:                                                            !
!    > L = real, range of x parameter space -L/2<= x <= +L/2           !
!    > N = integer, discretization of x parameter space: dim(H)=N+1    !
!    > omega = real, angular frequency                                 !
!    OUTPUT:                                                           !
!    > H = cmatrix, hamiltonian of the system                          !
!                                                                      !
! ---------------------------------------------------------------------!
! SUBROUTINES: none                                                    !
! ---------------------------------------------------------------------!
! ---------------------------------------------------------------------!
  use cmatrices

  contains

  function qho_H_init(L,N,omega) result(H)
    real          :: L,omega
    integer       :: N, ii

    real*16, dimension(:,:), allocatable :: elem_real
    type(cmatrix)                        :: H

    allocate(elem_real(N+1,N+1))
    elem_real = 0 * eleam_real

    ! diagonal
    do ii=1, N+1, 1
      elem_real(ii,ii) = ( 2  * (N*N)/(L*L) ) + omega*omega*((ii-1)*L/N - L/2)*((ii-1)*L/N - L/2)
    end do

    do ii=2, N+1, 1
      elem_real(ii,ii-1) = - (N*N)/(L*L)
      elem_real(ii-1,ii) = - (N*N)/(L*L)
    end do

    elem_real = 0.5* elem_real

    H = cmatrix_init(cmplx(X=elem_real,KIND=8))
  end function qho_H_init
end module qho

module debug
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
end module debug

program shroedingertimeindependent
  use cmatrices
  use qho

  real                       :: L, omega
  integer                    :: ii, N   ! number of discretized points
                                        ! -> dimension of H
  type(cmatrix)              :: H
  character(20)              :: folder
  integer                    :: iostat ! Checking types in READ(*,*)
  external::zgeev

  real*8, dimension(:), allocatable            :: eigenvalues
  complex(kind=8), dimension(:,:), allocatable :: eigenvectors
  double precision, dimension(:), allocatable  :: spacings

  integer :: infoeigens
  logical :: DEBUG

  DEBUG = .FALSE.

  print*, "--------------------------------------------"
  print*, "+        QUANTUM HARMONIC OSCILLATOR       +"
  print*, "--------------------------------------------"
  write(*,"(A)",advance='no') " + Type: L, N, omega and folder name: "
  read (*,*, iostat=iostat) L,N,omega,folder

  ! CHECK IF PARAMETERS ARE THE RIGHT TYPES
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!"
    print*, "+ N = int"
    print*, "+ L = real, int"
    print*, "+ ω = real, int"
    print*, "+ Exiting..."
    print*, "--------------------------------------------"
    stop
  end if

  ! CHECK IF PARAMETERS ARE IN THE RIGHT RANGE
  if(L <= 0 .OR. N<3 .OR. omega <=0) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!"
    print*, "+ (L > 0, N > 1, ω > 0)"
    print*, "+ Exiting..."
    print*, "--------------------------------------------"
    stop
  end if

  print*, "+ Data will be saved in: ./"//trim(folder)
  print*, "+ Lenght of x space (L): ", L
  print*, "+ Number of points  (N): ", N
  print*, "+ Angular frequency (ω): ", omega
  print*, "--------------------------------------------"
  print*, ""
  print*, "+ Computing the Hamiltonian..."

  H = qho_H_init(L,N,omega)

  if(N<3) then
    ! If we are in debug mode print also the imaginary part of H
    ! H should be real only
    if(DEBUG) then
      call cmatrix_print(H)
    else
      call cmatrix_print_real(H)
    end if
  else
    print*, " !!! H matrix is too big to be printed on screen !!!"
  end if

  print*, ""
  print*, "+ Computing Eigenvalues & Eigenvectors..."

  allocate(eigenvalues(N+1))
  allocate(eigenvectors(N+1,N+1))

  call cmatrix_herm_eigens(H,eigenvalues,eigenvectors,infoeigens)
  ! CHECK IF EIGENVALUE/EIGENVECTORS ARE COMPUTED PROPERLY
  if(DEBUG) then
    if(infoeigens == 0) then
      print*, "  [OK]"
    else
      print*, "  Something went wrong here"
    end if
  end if

  print*, "  Eigenvalues:"
  if(N<5) then
    do ii=1, N+1, 1
      if(DEBUG) then ! If we are in debug mode print also the imaginary part
                     ! that MUST BE 0
        print*, eigenvalues(ii)
      else
        print*, " ", real(eigenvalues(ii))
      end if
    end do
  else
    do ii=1, 5, 1
      if(DEBUG) then ! If we are in debug mode print also the imaginary part
                     ! that MUST BE 0
        print*, eigenvalues(ii)
      else
        print*, " ", real(eigenvalues(ii))
      end if
    end do
  end if

  call system("mkdir -p "//folder)
  open(42, file="./"//trim(folder)//"/eigenvalues.csv")
  print*, "  writing on file:"//" ./"//trim(folder)//"/eigenvalues.csv"
  do nn = 1, N+1, 1
    write(42, *) REAL(eigenvalues(nn))
  end do
  close(42)

  open(43, file="./"//trim(folder)//"/eigenvectors.csv")
  print*, "  writing on file:"//" ./"//trim(folder)//"/eigenvectors.csv"
  do nn = 1, N+1, 1
    write(43, *) REAL(eigenvectors(nn, :))
  end do
  close(43)

  deallocate(eigenvalues)
  deallocate(eigenvectors)

  print*, ""
  print*, "Done!"

end program shroedingertimeindependent
