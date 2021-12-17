module isingmod
  implicit none
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + ising_init_H:                                                      !
  !     Initialize the Hamiltonian of the ising model                    !
  !     INPUT:                                                           !
  !       > N = integer, number of spins                                 !
  !       > lambda = double precision, strenght of the external field    !
  !     OUTPUT:                                                          !
  !       > H = double precision, dimension(:,:), Hamiltonian matrix     !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  contains

  function ising_init_H(N,lambda) result(H)
    integer   :: N
    double precision :: lambda
    double precision, dimension(:,:), allocatable :: H, int_A, int_B

    integer :: ii,jj,kk

    allocate(H(2**N,2**N))
    H = 0.d0

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
      int_A = 0.d0
      int_B = 0.d0
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

      H = H - matmul(int_B,int_A)
      deallocate(int_A,int_B)
    end do
  end function
end module

module tensor_prod
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! FUNCTIONS:                                                           !
  ! + mat_tensor_I:                                                      !
  !     Performs tensor product between a matrix and an identity matrix  !
  !     of equal dimension                                               !
  !     INPUT:                                                           !
  !       > mat = double precision, dimension(:,:), matrix               !
  !     OUTPUT:                                                          !
  !       > matL = double precision, dimension(:,:), output matrix       !
  !                                                                      !
  ! + I_tensor_mat:                                                      !
  !     Performs tensor product between an identity matrix and a matrix  !
  !     of equal dimension                                               !
  !     INPUT:                                                           !
  !       > mat = double precision, dimension(:,:), matrix               !
  !     OUTPUT:                                                          !
  !       > matL = double precision, dimension(:,:), output matrix       !
  !                                                                      !
  ! + tens_prod:                                                         !
  !     Performs the tensor product between the two matrices             !
  !     INPUT:                                                           !
  !       > A, B = double precision, dimension(:,:), matrices            !
  !     OUTPUT:                                                          !
  !       > AoB = double precision, dimension(:,:), output matrix        !
  !                                                                      !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  contains

  function mat_tensor_I(mat) result(matL)
    double precision, dimension(:,:)              :: mat
    double precision, dimension(:,:), allocatable :: matL
    integer :: dimmat
    integer :: ii,jj ! coordinates of matrices
    integer :: kk

    if(size(mat,1)==size(mat,2)) then
      dimmat = size(mat,1)
      allocate(matL(dimmat**2,dimmat**2))
      matL = 0.d0

      do ii=0, dimmat-1, 1
        do jj=0, dimmat-1, 1
          do kk=1, dimmat, 1
            matL(dimmat*ii+kk,dimmat*jj+kk) = mat(ii+1,jj+1)
          end do
        end do
      end do
    end if
  end function

  function I_tensor_mat(mat) result(matR)
    double precision, dimension(:,:)              :: mat
    double precision, dimension(:,:), allocatable :: matR
    integer :: dimmat
    integer :: ii,jj ! coordinates of matrices
    integer :: kk

    if(size(mat,1)==size(mat,2)) then
      dimmat = size(mat,1)
      allocate(matR(dimmat**2,dimmat**2))
      matR = 0.d0

      do ii=1, dimmat, 1
        do jj=1, dimmat, 1
          do kk=0, dimmat-1, 1
            matR(ii+kk*dimmat,jj+kk*dimmat) = mat(ii,jj)
          end do
        end do
      end do
    end if
  end function

  function tens_prod(A, B) result(AoB)
    ! Computes A (X) B
    double precision, dimension(:,:) :: A, B
    double precision, dimension(:,:), allocatable :: AoB
    integer :: aa, bb, dimA, dimB

    dimA = size(A, 1)
    dimB = size(B, 1)
    allocate(AoB(dimA*dimB, dimA*dimB))

    do aa = 1, dimA, 1
      do bb = 1, dimB, 1
        AoB(dimB*(aa-1)+1:dimB*aa, dimB*(bb-1)+1:dimB*bb ) = A(aa, bb)*B
      end do
    end do
  end function tens_prod
end module tensor_prod

module rg
  use tensor_prod
  contains
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  ! SUBROUTINES:                                                         !
  ! + init_interaction_H                                                 !
  !     Initialize the interaction matrix between the two subsystems     !                                                    !
  !       > N: number of spins for each subsystems                       !
  !       > HL, HR: double precision, dimension(:,:), interaction matrix !
  !                                                                      !
  ! + diagonalize_H                                                      !
  !     Diagonalize the hamiltonian and get the projections matrix       !
  !       > H: double precision, dimension(:,:), hamiltonian of size 2N  !
  !       > eigenvalues: double precision, dimension(:), eigenvalues     !
  !       > nvec: integer, number of eigenvector to store                !
  !       > eigenvectors: double precision, dimension(:,:), projection   !
  ! ---------------------------------------------------------------------!
  ! ---------------------------------------------------------------------!
  subroutine init_interaction_H(N,HL,HR)
    double precision, dimension(:,:), allocatable :: HL, HR
    integer :: kk, N

    allocate(HL(2**(N),2**(N)),HR(2**(N),2**(N)))
    HL = 0.d0
    HR = 0.d0

    do kk = 0, 2**(N-1)-1, 1
      HL(2+(2*kk),1+(2*kk)) = 1
      HL(1+(2*kk),2+(2*kk)) = 1

      HR(2**(N-1)+1+kk,1+kk) = 1
      HR(1+kk,2**(N-1)+1+kk) = 1
    end do
  end subroutine

  subroutine diagonalize_H(H, eigenvalues, nvec, eigenvectors)
    double precision, dimension(:,:), allocatable :: H
    double precision, dimension(:), allocatable   :: eigenvalues
    double precision, dimension(:,:), allocatable :: eigenvectors

    double precision, dimension(:,:), allocatable              :: supp
    integer :: LWORK, LIWORK, INFO, N, M, LWMAX, dims(2), nvec
    double precision :: VL, VU, ABSTOL
    double precision, dimension(:), allocatable                :: WORK
    integer, dimension(:), allocatable                         :: ISUPPZ, IWORK

    dims = shape(H)
    N    = dims(1)
    LWMAX = 100000
    ABSTOL = 0.
    VU = 0.
    VL = 0.

    allocate(WORK(LWMAX))
    allocate(IWORK(LWMAX))
    allocate(ISUPPZ(2*nvec))
    allocate(eigenvectors(N, nvec))
    allocate(supp(N, N))
    allocate(eigenvalues(N))

    supp = H
    LWORK = -1
    LIWORK = -1

    ! Computing optimal workspace...
    call DSYEVR('V', 'I', 'U', N, supp, N, VL, VU, 1, nvec, &
                ABSTOL, M, eigenvalues, eigenvectors, N, ISUPPZ, &
                WORK, LWORK, IWORK, LIWORK, INFO)
    LWORK = min(LWMAX, int(WORK(1)))
    LIWORK = min(LWMAX, int(IWORK(1)))

    ! Computing eigenvalues...
    call DSYEVR('V', 'I', 'U', N, supp, N, VL, VU, 1, nvec, &
                ABSTOL, M, eigenvalues, eigenvectors, N, ISUPPZ, &
                WORK, LWORK, IWORK, LIWORK, INFO)

    deallocate(IWORK, WORK, ISUPPZ)
  end subroutine diagonalize_H
end module rg

module other
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
  ! ---------------------------------------------------------------------!
  ! TYPES: none                                                          !
  ! ---------------------------------------------------------------------!
  contains

  character(len=20) function str(k)
    ! "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  subroutine project(P, M, Mred)
    double precision, dimension(:,:)              :: P, M
    double precision, allocatable, dimension(:,:) :: Mred

    if(allocated(Mred)) then
      deallocate(Mred)
    end if

    allocate( Mred( size(P,2), size(P,2) ) )
    Mred = matmul( matmul( transpose(P), M), P )
  end subroutine project

  function identity(dim) result(id)
    integer                          :: dim, ii
    double precision, dimension(:, :), allocatable :: id

    allocate(id(dim, dim))
    id = 0d0

    do ii = 1, dim
      id(ii, ii) = 1d0
    end do
    end function identity
end module other

program ising_rg
  use other
  use tensor_prod
  use isingmod
  use rg

  integer          :: N, nit ! initial number of spins, max number of iterations
  double precision :: lambda ! strenght of external field

  integer :: it ! iteration of the algorithm it goes from 1 to nit

  ! DEBUGGING PURPOSES
  logical :: DEBUG
  integer :: iostat

  ! of each iteration...
  double precision, dimension(:,:), allocatable :: HN         ! first hamiltonian
  double precision, dimension(:,:), allocatable :: H2N        ! doubled hamiltonian
  double precision, dimension(:,:), allocatable :: HL, HR     ! interac hamiltonian
  double precision, dimension(:,:), allocatable :: HLred, HRred     ! interac hamiltonian
  double precision, dimension(:,:), allocatable :: P          ! Project

  double precision, dimension(:), allocatable :: evls        ! eigenvalues
  double precision                            :: oldevl      ! eigenvalue of the previous iteration
  logical                                     :: converged   ! check if already converged
  integer*8                                   :: sizeofspace ! number of spins at each iterations
  integer                                     :: stopprog

  character(20) :: folder, file


  DEBUG = .TRUE.
  converged = .FALSE.

  print*, "+------------------------------------------+"
  print*, "+   renormalization group algorith on the  +"
  print*, "+                ISING MODEL               +"
  print*, "+------------------------------------------+"

  print*, "+ This program applies the rg-algorithm to +"
  print*, "+ the Ising Model                          +"
  print*, "+------------------------------------------+"

  print*, "+ N: starting Number of spin1/2 particles  +"
  print*, "+ nit: number of iterations                +"
  print*, "+ lambda: intensity of external field      +"
  print*, "+ stopprog: stop the program if convergence+"
  print*, "+           has been reached               +"
  print*, "+           stop -> 1                      +"
  print*, "+           dont -> 0                      +"
  print*, "+------------------------------------------+"

  write(*,"(A)",advance='no') " + Type N, nit, lambda, stopprog, folder, file: "
  read (*,*, iostat=iostat) N, nit, lambda, stopprog, folder, file

  call system("mkdir -p ./data/"//trim(folder))

  ! Store the gs energy at each iteration
  open(1, file = "./data/"//trim(folder)//"/"//trim(file)//".csv")

  ! Store the number of iterations to each stability
  open(2, file = "./data/"//trim(folder)//"/"//trim(file)//"convergence.csv")

  ! Check if parameters are the right types
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!          +"
    print*, "+ N = integer                              +"
    print*, "+ nit = integer                            +"
    print*, "+ lambda = float                           +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check if parameters are in the right ranges
  if(N < 2 .OR. nit < 1) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!         +"
    print*, "+ (N > 1)                                  +"
    print*, "+ (nit > 0)                                +"
    print*, "+ Exiting...                               +"
    print*, "+------------------------------------------+"
    print*, "+------------------------------------------+"
    stop
  end if

  ! Check on READ
  if(DEBUG) then
    write(*,"(A)",advance='no') " + N:     "
    print*, N
    write(*,"(A)",advance='no') " + nit:   "
    print*, nit
    write(*,"(A)",advance='no') " + lambda:"
    print*, lambda
  end if
  print*, "+                                          +"

  HN = ising_init_H(N,lambda)
  sizeofspace = N

  if(DEBUG) then
    print*, "+ First Hamiltonian initialized            +"
  end if

  call init_interaction_H(N,HL,HR)
  print*, "+ Iterating...                             +"
  ! iterations
  do it = 1, nit, 1
    if(modulo(it*10,nit)==0) then
      write(*,'(A,I6,A,I6,A)',advance='yes') " + ", it, "  /", nit, "                          +"
    end if

    allocate(H2N(2**(2*N),2**(2*N)))

    sizeofspace = 2*sizeofspace

    H2N = mat_tensor_I(HN) + I_tensor_mat(HN) + tens_prod(HL,HR)
    HLred = tens_prod(HL, identity(2**N))
    HRred = tens_prod(identity(2**N), HR)

    call diagonalize_H(H2N, evls, 2**N, P)

    call project(P, H2N, HN)
    call project(P, HLred, HL)
    call project(P, HRred, HR)

    if(abs(oldevl - evls(1)/sizeofspace)<1d-12) then
      if(converged .eqv. .FALSE.) then
        print*, "+ Algorithm has reached convergence"
        if(stopprog == 1) then
          converged = .TRUE.
          write(2,*) it
          stop
        end if
      end if
    end if

    write(1,*) evls(1)/sizeofspace
    oldevl = evls(1)/sizeofspace

    deallocate(H2N,evls,P)

  end do
  if(converged .eqv. .FALSE.) then
    print*, " Algorithm has NOT reached convergence     +"
    write(2,*) -1
  end if
  print*, "+------------------------------------------+"
  print*, "+ Done!                                    +"
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program ising_rg
