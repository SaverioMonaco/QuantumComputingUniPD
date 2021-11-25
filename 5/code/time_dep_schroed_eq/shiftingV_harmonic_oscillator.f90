! Quantum Information and Computing 2021-2022
! Saverio Monaco, MAT. 2012264
! Week 5: Solving time dependent Schroedinger Equation: Shifting V QHO

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
    cmat%element = cmplx(elem_real,elem_imag,kind=8)
    cmat%trace   = cmatrix_trace(cmat)
    cmat%det     = 0
  end function cmatrix_randinit

  function cmatrix_randinit_hermitian(n,range) result(cmat)
    integer                              :: n
    real                                 :: range
    real*16, dimension(:,:), allocatable :: elem_real,elem_imag
    integer                              :: ii,jj
    type(cmatrix)                        :: cmat

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
    cmat%element = cmplx(elem_real,elem_imag,kind=8)
    cmat%trace   = cmatrix_trace(cmat)
    cmat%det     = 0
  end function cmatrix_randinit_hermitian

  function cmatrix_init(array2d) result(cmat)
    complex*16, dimension(:,:)  :: array2d
    type(cmatrix)               :: cmat

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
    elem_real = 0

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

module qho_split_op
  use qho
  use cmatrices

  contains
  subroutine qho_split_U1(psix0, psix1, current_t, T, N, Lx)
    complex(kind=8), dimension(:)              :: psix0, psix1
    real                                       :: current_T, T, Lx
    real*8, dimension(N)                       :: x
    integer                                    :: ii, N

    do ii=0,N,1
      x(ii+1) = -0.5d0*Lx + Lx/N*ii
    end do

    psix1 = psix0*( EXP(dcmplx(0d0,-0.5d0)*(T/N)*(.5*(x-current_t/T)**2) ) )
  end subroutine qho_split_U1

  subroutine qho_split_fourier(psi_in, psi_out)
    complex(kind=8), dimension(:)     :: psi_in, psi_out
    ! FFTW Related
    integer*8                         :: plan

    call dfftw_plan_dft_1d(plan, size(psi_in), psi_in, psi_out, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute(plan)!, signal, signal_out)

  end subroutine qho_split_fourier

  subroutine qho_split_U2(psik1, psik2, T, N, Lp)
    complex(kind=8), dimension(:)              :: psik1, psik2
    real                                       :: T, Lp
    real*8, dimension(N)                       :: p
    integer                                    :: ii, N

    do ii=0,N,1
      p(ii+1) = -0.5d0*Lp + Lp/N*ii
    end do

    psik2 = psik1*( EXP(dcmplx(0d0,-0.5d0)*(T/N)*((p)**2) ) )
  end subroutine qho_split_U2

end module qho_split_op

module debugmod
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

end module debugmod

program shroedingertimedependent
  use debugmod
  use cmatrices
  use qho
  use qho_split_op

  external::zheev

  type(cmatrix)              :: H ! Initial Hamiltonian for t=0
  real                       :: Lx, Lp, T, current_t
  integer                    :: lev, N, ii, tt
  integer                    :: Nt      ! Discretization of time
  character(20)              :: folder
  integer                    :: iostat  ! Checking types in READ(*,*)

  real*8, dimension(:), allocatable            :: eigenvalues
  complex(kind=8), dimension(:,:), allocatable :: eigenvectors
  ! Just the eigenvector of the energy level we are interested in at t=0
  complex(kind=8), dimension(:,:), allocatable :: eigenvector

  ! Vectors for:
  ! tmp_k1: after Ux, Fourier transformed
  ! tmp_k2: after Uk
  ! tmp_x2: after inverse transforming
  complex(kind=8), dimension(:), allocatable   :: tmp_x1, tmp_k1, tmp_k2, tmp_x2, tmp_x3

  integer :: infoeigens
  logical :: DEBUG

  DEBUG = .FALSE.

  print*, "--------------------------------------------"
  print*, "+   TIME DEPENDENT SCHROEDINGER EQUATION   +"
  print*, "+------------------------------------------+"
  print*, "+ Lx:  Lenght of x space                   +"
  print*, "+ Lp:  Lenght of p space                   +"
  print*, "+ N:   Number of points                    +"
  print*, "+ T:   Total time                          +"
  print*, "+ Nt:  Number of time-points               +"
  print*, "+ lev: Energy level to evolve              +"
  print*, "+------------------------------------------+"
  print*, "+ Type: Lx, Lp, N, T, Nt, lev and folder:  +"
  write(*,"(A)",advance='no') " + "

  read (*,*, iostat=iostat) Lx, Lp, N, T, Nt, lev, folder

  ! CHECK IF PARAMETERS ARE THE RIGHT TYPES
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!"
    print*, "+ N   = int"
    print*, "+ Nt  = int"
    print*, "+ Lx  = real, int"
    print*, "+ Lp  = real, int"
    print*, "+ T   = real, int"
    print*, "+ lev = int"
    print*, "+ Exiting..."
    print*, "--------------------------------------------"
    stop
  end if

  ! CHECK IF PARAMETERS ARE IN THE RIGHT RANGE
  if(Lx <= 0 .OR. Lp <= 0 .OR. T<=0 .OR. Nt<2 .OR. N<3) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!"
    print*, "+ (Lx > 0, Lp > 0, N > 2, Nt > 1, T > 0, lev >= 0)"
    print*, "+ Exiting..."
    print*, "--------------------------------------------"
    stop
  end if

  print*, "+ Data will be saved in: ./"//trim(folder)
  print*, "+ Lenght of x space     (Lx): ", Lx
  print*, "+ Lenght of p space     (Lp): ", Lp
  print*, "+ Number of points       (N): ", N
  print*, "+ Total time             (T): ", T
  print*, "+ Number of time-points (Nt): ", Nt
  print*, "+ Energy level         (lev): ", lev
  print*, "+------------------------------------------+"

  print*, "+"
  print*, "+ Computing eigenstate for t = 0..."
  write(*,"(A)",advance='no') " +  ∟ Computing the Hamiltonian "
  H = qho_H_init(Lx,N,1.0)
  print*, "           [DONE]"

  if(N<4) then
    ! If we are in debug mode print also the imaginary part of H
    ! H should be real only
    if(DEBUG) then
      print*, "+ "
      call cmatrix_print(H)
    end if
  end if

  write(*,"(A)",advance='no') " +  ∟ Computing Eigenvalues & Eigenvectors "
  allocate(eigenvalues(N+1))
  allocate(eigenvectors(N+1,N+1))
  call cmatrix_herm_eigens(H,eigenvalues,eigenvectors,infoeigens)
  allocate(eigenvector(N+1, Nt+1))
  eigenvector(:,1) = eigenvectors(:,lev+1)

  print*, "[DONE]"

  if(DEBUG) then
    print*, "+"
    print*, "+ Eigenvalues:"
    do ii=1,5,1
      print*, "+ ", eigenvalues(ii)
    end do
    print*, "+    ..."
  end if


  if(DEBUG) then
    print*, "+"
    print*, "+ Writing initial state wavefunction at:"
    print*, "+   ./"//trim(folder)//"/DEBUG_PSI0.csv"

    call system("mkdir -p "//folder)

    open(42, file="./"//trim(folder)//"/DEBUG_PSI0.csv")
    do ii=1,N+1,1
      write(42, *) eigenvector(ii,1)
    end do
    close(42)
    print*, "+ [DONE]"
  end if

  ! Free some memory
  !deallocate(eigenvectors)

  if(DEBUG) then
    print*, "+"
    print*, "+ Psi(x) lev:", lev
    do ii=1,5,1
      print*, "+ ", eigenvectors(ii,lev+1)
    end do
    print*, "+   ..."
  end if

  ! ###########################################################################
  ! ########################## SPLIT OPERATOR METHOD ##########################
  ! ###########################################################################
  print*, "+"
  print*, "+ Evolving:"
  allocate(tmp_x1(N+1), tmp_k1(N+1), tmp_k2(N+1), tmp_x2(N+1), tmp_x3(N+1) )
  current_t = 0

  do tt=1, Nt+1, 1
    if(mod(tt,Nt/10)==0) then
      write(*,*) '+ ', tt*100/Nt, '%'
    end if

    current_t = current_t + T/Nt

    ! Evolve Ux Ψx
    call qho_split_U1(eigenvector(:,tt), tmp_x1, current_t, T, N, Lx)

    ! Transform Ψx -> Ψk
    call qho_split_fourier(tmp_x1, tmp_k1)

    ! Evolve Uk Ψk
    call qho_split_U2(tmp_k1, tmp_k2, T, N, Lp)

    ! Transform Ψx -> Ψk
    call qho_split_fourier(tmp_k2, tmp_x2)

    ! Evolve Ux Ψx
    call qho_split_U1(tmp_x2, tmp_x3, current_t, T, N, Lx)

    eigenvector(:,tt+1) = tmp_x3
  end do
  print*, "+ [DONE]"
  ! ###########################################################################
  ! ###################### END OF SPLIT OPERATOR METHOD #######################
  ! ###########################################################################

  print*, "+"
  print*, "+------------------------------------------+"
  print*, "+ Writing on files: ./"//trim(folder)//"/*_wavefunc.csv"

  call system("mkdir -p "//folder)

  open(22, file="./"//trim(folder)//"/REAL_wavefunc.csv")
  open(23, file="./"//trim(folder)//"/IMAG_wavefunc.csv")
  do tt=1, Nt+1, 1
    write(22, *) REAL(eigenvector(:,tt))
    write(23, *) AIMAG(eigenvector(:,tt))
  end do
  close(22)
  close(23)

  print*, "+ [DONE]"
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program shroedingertimedependent

! Compile:
!  gfortran shiftingV_harmonic_oscillator.f90 -o sh_qho -llapack -Wall -lfftw3 -frecursive
