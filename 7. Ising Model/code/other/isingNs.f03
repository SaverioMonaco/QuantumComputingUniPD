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
  use isingmod

  implicit none

  integer          :: N, d, SEP, Nmax
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
  DEBUG = .FALSE.

  print*, "--------------------------------------------"
  print*, "+                ISING MODEL               +"
  print*, "+------------------------------------------+"
  print*, "+ This program investigates quantities for +"
  print*, "+ the Ising Model                          +"
  print*, "+------------------------------------------+"
  print*, "+ Nmax: Number of spin 1/2 particles          +"
  print*, "+ λ: Number of lambdas (range is [0,3])+"
  print*, "+------------------------------------------+"
  write(*,"(A)",advance='no') " + Type N, λ, folder: "
  read (*,*, iostat=iostat) N, lambda, folder

  call system("mkdir -p ./data/Ns/"//trim(folder))
  open(1, file = "./data/Ns/"//trim(folder)//'/cputime')
  open(2, file = "./data/Ns/"//trim(folder)//'/eigenvalue')

  ! Check if parameters are the right types
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!          +"
    print*, "+ Nmax = integer                           +"
    print*, "+ λ    = float                             +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check if parameters are in the right ranges
  if(N < 2) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!         +"
    print*, "+ (Nmax > 1)                               +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check on READ
  if(DEBUG) then
    write(*,"(A)",advance='no') " + Nmax:"
    print*, Nmax
    write(*,"(A)",advance='no') " + d:"
    print*, d
    write(*,"(A)",advance='no') " + λ:"
    print*, lambda
    if(SEP==1) then
      print*, "+ System is separable"
    else
      print*, "+ System is NOT separable"
    end if
  end if

  do N = 2, Nmax, 1
    ! Preparation

    ! Building the Hamiltonian
    call cpu_time(cpu_start)  ! start time

    ! Diagonalizing and stuff
    H = ising_init_H(N,lambda)
    if(DEBUG) then
      print*, "+ N:", lambda
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

    if(DEBUG) then
      write(*,'(A)',advance='no') " +   "
      do ii = 1, 3, 1
        write(*,'(A,ES0.2,A)', advance='no') "  ", eigenv(d**N - ii + 1), ", "
      end do
      print*, ""
    end if

    deallocate(RWORK, VR, H)

    call cpu_time(cpu_finish) ! end time

    write(1,*) cpu_finish - cpu_start
    write(2,*) eigenv

    deallocate(eigenv)
  end do

  print*, "+ Saving results in:                       +"
  print*, "+   ./data/Ns/"//trim(folder)
  print*, "+------------------------------------------+"


  print*, "+ All done!                                +"
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program ising
