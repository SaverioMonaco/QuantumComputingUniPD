module basechange
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + basechange_to                                                      !
  !   INPUT:                                                             !
  !   > b_to = integer, base for the conversion 10 -> b_to               !
  !   > number = integer, number to convert number_10 -> number_b_to     !
  !   > N = integer, size of number: 4_2 is 0100 if N = 4                !
  !   OUTPUT:                                                            !
  !   > number_b_to = integer(N), number converted                       !
  !                                                                      !
  ! + basechange_from                                                    !
  !   INPUT:                                                             !
  !   > b_from = integer, base for the conversion b_to -> 10             !
  !   > number_from = integer(N), number to convert from b_from to 10    !
  !   > N = integer, size of number: 4_2 is 0100 if N = 4                !
  !   OUTPUT:                                                            !
  !   > number_b10 = integer, number converted                           !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES: none                                                    !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  contains

  function basechange_to(b_to, number, N) result(number_b_to)
    integer :: b_to, number, N, ii
    integer, dimension(N) :: number_b_to

    number_b_to = 0*number_b_to ! Allocate to 0, just to be sure
    do ii = 1, N, 1
      number_b_to(N - ii + 1) = modulo(number, b_to)
      number = number/b_to
    end do
  end function

  function basechange_from(b_from, number_from, N) result(number_b10)
    integer :: b_from, number_b10, N, ii
    integer, dimension(N) :: number_from

    do ii = 1, N, 1
      number_b10 = number_b10 + number_from(N - ii + 1)*b_from**(ii - 1)
    end do
  end function
end module

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
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(2))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(4))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(8))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (real(4))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (real(8))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (logical)
            print*, "+ ", msg, " => [OK], Variable = ", content
        end select
      else
        print*, "+ ", msg, " => [OK]"
      end if
    end if
    end subroutine
end module debugmod

module densmat
  use basechange
  implicit none
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES:                                                               !
  ! + qsystem:                                                           !
  !           Structure for quantum multiple subsystems                  !
  !           > N           : integer, number of states                  !
  !           > d           : actual, dimension of each Hilbert space    !
  !           > separability: wether the system is separable or not      !
  !           > waves       : elements/coefficients of the systems       !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + densmat_pure_init                                                  !
  !   Initializing function for qsystems                                 !
  !   INPUT:                                                             !
  !   > N: integer, number of states                                     !
  !   > d: integer, dimension of each Hilbert spaces                     !
  !   > SEP: integer if separable (SEP=1), if not (SEP=0)                !
  !   > DEBUG_: logical, prints bytes allocated (can be too many)        !
  !   OUTPUT:                                                            !
  !   > PHI: qsystem                                                     !
  !                                                                      !
  ! + densmat_computerho1                                                !
  !   Compute partial trace                                              !
  !   INPUT                                                              !
  !   > rho: double complex, dimension(d**2,d**2), density matrix        !
  !   > d: integer, size of each Hilbert space                           !
  !   OUTPUT:                                                            !
  !   > rho1: double complex, dimension(d,d), density matrix             !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES:                                                         !
  ! + densmat_genstates                                                  !
  !   Generate random numbers (dxN or d^N double complex numbers)        !
  !   INPUT:                                                             !
  !   > PHI: qsystem                                                     !
  !                                                                      !
  ! + densmat_readcoeffs                                                 !
  !   Read coefficient from terminal input                               !
  !   INPUT:                                                             !
  !   > PHI: qsystem                                                     !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  type qsystem
    integer :: N ! Number of systems
    integer :: d ! number of states

    logical :: separability ! Wether the whole system is separable
                            ! (no interactions) or not, the number of
                            ! coefficients changes:
                            ! separability == TRUE  -> #C = d X N
                            !              == FALSE -> #C = d ^ N

    double complex, dimension(:), allocatable   :: waves
    double complex, dimension(:,:), allocatable :: rho
  end type qsystem

  contains

  function densmat_pure_separable(N, d) result(PHI)
    integer   :: N, d
    type(qsystem) :: PHI

    PHI%N = N
    PHI%d = d
    PHI%separability = .TRUE.

    allocate(PHI%waves(N*d))
  end function

  function densmat_pure_unseparable(N, d) result(PHI)
    integer   :: N, d
    type(qsystem) :: PHI

    PHI%N = N
    PHI%d = d
    PHI%separability = .FALSE.

    allocate(PHI%waves(d**N))
  end function

  function densmat_pure_init(N, d, SEP, DEBUG) result(PHI)
    integer       :: N, d, SEP
    integer*8     :: nbytes
    type(qsystem) :: PHI
    logical       :: DEBUG

    if(SEP == 1) then
      PHI = densmat_pure_separable(N, d)
      nbytes = 16*(N*d)
      if(DEBUG) then
        print*, "+ Allocating", nbytes, "bytes"
      end if
    else if(SEP == 0) then
      PHI = densmat_pure_unseparable(N, d)
      nbytes = 16*(d**N)
      if(DEBUG) then
        print*, "+ Allocating", nbytes, "bytes"
      end if
    else
      print*, "+ Unvalid state, exiting .... "
      stop
    end if
  end function

  subroutine densmat_genstates(PHI)
    type(qsystem) :: PHI
    double precision, dimension(:), allocatable :: psi_real, psi_imag
    double precision                            :: norm

    if(PHI%separability .eqv. .TRUE.) then
      ! If separable allocate dXN
      allocate(psi_real(PHI%d*PHI%N),psi_imag(PHI%d*PHI%N))
    else
      ! else allocate d^N
      allocate(psi_real(PHI%d**PHI%N),psi_imag(PHI%d**PHI%N))
    end if

    call random_number(psi_real) ! Generate all real
    call random_number(psi_imag) ! Generate all imag
    PHI%waves = dcmplx(psi_real, psi_imag)
    ! Normalizing
    norm = SUM(ABS(PHI%waves(:))**2)
    ! Assigning
    PHI%waves(:) =  PHI%waves(:)/(sqrt(norm))
  end subroutine

  subroutine densmat_readcoeffs(PHI)
    type(qsystem) :: PHI
    double precision :: real_coeff, imag_coeff, norm
    integer :: ii,jj

    print*, "+ Assign coefficiens:"
    do jj = 0, PHI%d**PHI%N - 1, 1
      ii = jj
      write(*,"(A)",advance='no') " + Coefficient for state"
      write(*,"(A)",advance='no') "  |"
      write(*,'(*(I4))', advance='no') basechange_to(PHI%d,ii,PHI%N)
      write(*,"(A)",advance='no') " > : (real imag)  "
      read (*,*) real_coeff, imag_coeff
      PHI%waves(jj + 1) = dcmplx(real_coeff,imag_coeff)
    end do

    write(*,"(A)",advance='no') " + Normalizing..."
    ! Normalizing
    norm = SUM(ABS(PHI%waves(:))**2)
    ! Assigning
    PHI%waves(:) =  PHI%waves(:)/(sqrt(norm))
    print*, "  Done!"

  end subroutine

  function densmat_computerho1(rho,d) result(rho1)
    integer :: d
    double complex, dimension(:,:) :: rho
    double complex, dimension(:,:), allocatable :: rho1
    integer :: ii,jj,kk

    allocate(rho1(d,d))

    do ii = 1, d
      do jj = 1, d
        rho1(ii,jj) = 0
        do kk = 0, d - 1
          rho1(ii,jj) = rho1(ii,jj) + rho(d*(ii-1) + 1 + kk, d*(jj-1) + 1 + kk)
        end do
      end do
    end do
  end function
end module densmat

module isingmod
  implicit none
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  contains

  function ising_init_H(N,lambda) result(H)
    integer   :: N
    double precision :: lambda
    double complex, dimension(:,:), allocatable :: H, int_A, int_B

    integer :: ii,jj,kk,ll

    allocate(H(2**N,2**N))
    H = 0.0 * H

    ! External field part: \lambda \sum_i^N \sigma_z^i
    do ii = 1, N, 1
      print*, "stop"
      do jj = 1, 2**N, 1
        H(jj,jj) = H(jj,jj) + -2*(modulo( (jj-1)/int(2**(N-ii)),2) ) +1
      end do
    end do
    H = lambda * H ! Adding the magnetization field factor

    ! Interaction part -\sum_i^{N-1}\sigma_x^{i+1}\sigma_x^i
    do ii = 1, N-1, 1
      allocate(int_A(2**N,2**N))
      allocate(int_B(2**N,2**N))
      int_A = int_A * 0.0
      int_B = int_B * 0.0
      do kk = 0,2**(ii-1)-1,1
        do jj=1,2**(N-ii),1
          int_A(kk*(2**(N-ii+1))  + 2**(N-ii)+jj, kk*(2**(N-ii+1))  + jj) = 1
          int_A(kk*(2**(N-ii+1))  + jj, kk*(2**(N-ii+1))  + 2**(N-ii)+jj) = 1
        end do
      end do
      do kk = 0,2**(ii)-1,1
        do jj=1,2**(N-ii-1),1
          int_B(kk*(2**(N-ii)) + 2**(N-ii-1)+jj, kk*(2**(N-ii)) + jj) = 1
          int_B(kk*(2**(N-ii)) + jj, kk*(2**(N-ii)) + 2**(N-ii-1)+jj) = 1
        end do
      end do
      if(.False. .eqv. .True.) then
      print*, "mata"
      do jj = 1, ubound(int_A, 1)
        print*, "|", real(int_A(jj, :)), "|"
      end do
      print*, "matB"
      do jj = 1, ubound(int_B, 1)
        print*, "|", real(int_B(jj, :)), "|"
      end do
      end if
      H = H - matmul(int_B,int_A)
      deallocate(int_A,int_B)
    end do
  end function
end module

program ising

  use debugmod
  use densmat
  use isingmod

  implicit none

  integer          :: N, d, SEP, nlams,ll
  double precision :: lambda
  character(20)    :: folder
  double complex, dimension(:,:), allocatable :: H

  ! Debugging
  logical :: DEBUG
  integer :: iostat,ii

  ! CPU times
  real*8 :: cpu_start, cpu_finish ! for the CPU times

  ! LAPACK variables
  double precision, dimension(:), allocatable   :: RWORK
  integer                                       :: INFO, LWORK
  integer, parameter                            :: LWMAX = 100000
  complex*16                                    :: WORK(LWMAX)
  complex(kind=8), dimension(:,:), allocatable  :: VR
  double precision, dimension(:), allocatable   :: eigenv

  ! PARAMETERS
  d   = 2 ! Spin 1/2 paticles
  SEP = 0
  DEBUG = .TRUE.

  print*, "--------------------------------------------"
  print*, "+                ISING MODEL               +"
  print*, "+------------------------------------------+"
  print*, "+ This program investigates quantities for +"
  print*, "+ the Ising Model                          +"
  print*, "+------------------------------------------+"
  print*, "+ N: Number of spin 1/2 particles          +"
  print*, "+ nlams: Number of lambdas (range is [0,3])+"
  print*, "+------------------------------------------+"
  write(*,"(A)",advance='no') " + Type N, nlams, folder: "
  read (*,*, iostat=iostat) N, nlams, folder

  call system("mkdir -p ./data/"//trim(folder))
  open(1, file = "./data/"//trim(folder)//'/cputime')
  open(2, file = "./data/"//trim(folder)//'/eigenvalue')

  ! Check if parameters are the right types
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!          +"
    print*, "+ N = integer                              +"
    print*, "+ nlams = integer                          +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check if parameters are in the right ranges
  if(N < 2) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!         +"
    print*, "+ (N > 1)                                  +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check on READ
  if(DEBUG) then
    write(*,"(A)",advance='no') " + N:"
    print*, N
    write(*,"(A)",advance='no') " + d:"
    print*, d
    write(*,"(A)",advance='no') " + nlams:"
    print*, nlams
    if(SEP==1) then
      print*, "+ System is separable"
    else
      print*, "+ System is NOT separable"
    end if
  end if

  do ll = 0, nlams, 1
    lambda = (ll)*3/float(nlams)
    ! Preparation

    ! Building the Hamiltonian
    call cpu_time(cpu_start)  ! start time

    ! Diagonalizing and stuff
    H = ising_init_H(N,lambda)
    if(DEBUG) then
      print*, "+ lambda:", lambda
      do ii = 1, ubound(H, 1)
        print*, "|", real(H(ii, :)), "|"
      end do
    end if

    allocate(RWORK(3*(d**N)-2))
    allocate(VR(d**N,d**N))
    allocate(eigenv(d**N))
    ! Compute optimal size of workspace
    LWORK = -1

    call ZHEEV('N', 'U', d**N, H, d**N, eigenv, WORK,LWORK,RWORK,INFO)
    LWORK = min(LWMAX, int(WORK(1)))

    ! Compute eigenvalues
    call ZHEEV('N', 'U', d**N, H, d**N, eigenv, WORK,LWORK,RWORK,INFO)

    write(*,'(A)',advance='no') " +   "
    do ii = 1, 3, 1
      write(*,'(A,ES0.2,A)', advance='no') "  ", eigenv(d**N - ii + 1), ", "
    end do
    print*, ""

    deallocate(RWORK, VR, H)

    call cpu_time(cpu_finish) ! end time

    write(1,*) cpu_finish - cpu_start
    write(2,*) eigenv

    deallocate(eigenv)
  end do

  print*, "+ Saving results in:                       +"
  print*, "+   ./data/"//trim(folder)
  print*, "+------------------------------------------+"


  print*, "+ All done!                                +"
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program ising
