module isingmod
  implicit none
  ! ---------------------------------------------------------------------!
  ! --------------------------- DOCUMENTATION ---------------------------!
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

  ! Performs I \ocircle sigma^x_N \ocircle sigma^x_{N+1} \ocircle I
  function init_interaction_H(N) result(Hin)
    double precision, dimension(:,:), allocatable :: Hin, HL, HR
    integer :: kk

    allocate(HL(2**(N),2**(N)),HR(2**(N),2**(N)))
    HL = 0.d0
    HR = 0.d0

    do kk = 0, 2**(N-1)-1, 1
      HL(2+(2*kk),1+(2*kk)) = 1
      HL(1+(2*kk),2+(2*kk)) = 1

      HR(2**(N-1)+1+kk,1+kk) = 1
      HR(1+kk,2**(N-1)+1+kk) = 1
    end do

    allocate(Hin(2**(2*N),2**(2*N)))

    Hin = tens_prod(HL,HR)
  end function

  subroutine findkEigenvalue(matrix, eigenvalues, kk, Z)

        double precision, dimension(:,:), allocatable, intent(IN)    :: matrix
        double precision, dimension(:), allocatable, intent(INOUT) :: eigenvalues
        double precision, dimension(:,:), allocatable, intent(OUT) :: Z

        double precision, dimension(:,:), allocatable              :: supp
        integer :: LWORK, LIWORK, INFO, N, M, LWMAX, dims(2), kk
        double precision :: VL, VU, ABSTOL
        double precision, dimension(:), allocatable                :: WORK
        integer, dimension(:), allocatable                         :: ISUPPZ, IWORK


        dims = shape(matrix)
        N    = dims(1)
        LWMAX = 100000
        ABSTOL = 0.
        VU = 0.
        VL = 0.

        allocate(WORK(LWMAX))
        allocate(IWORK(LWMAX))
        allocate(ISUPPZ(2*kk))
        allocate(Z(N, kk))
        allocate(supp(N, N))
        allocate(eigenvalues(N))

        supp = matrix

    !   compute optimal size of workspace
        LWORK = -1
        LIWORK = -1
        !print*, "Computing optimal workspace..."
        call DSYEVR('V', 'I', 'U', N, supp, N, VL, VU, 1, kk, &
                    ABSTOL, M, eigenvalues, Z, N, ISUPPZ, &
                    WORK, LWORK, IWORK, LIWORK, INFO)
        LWORK = min(LWMAX, int(WORK(1)))
        LIWORK = min(LWMAX, int(IWORK(1)))

    !   compute eigenvalues
        !print*, "Computing eigenvalues..."
        call DSYEVR('V', 'I', 'U', N, supp, N, VL, VU, 1, kk, &
                    ABSTOL, M, eigenvalues, Z, N, ISUPPZ, &
                    WORK, LWORK, IWORK, LIWORK, INFO)

        if (INFO .NE. 0) then
            print*, 'Failure!'
            stop
        end if

        deallocate(IWORK, WORK, ISUPPZ)

    end subroutine findkEigenvalue
end module rg

module other
  contains

  character(len=20) function str(k)
    ! "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
end module other

program ising_rg
  use other
  use tensor_prod
  use isingmod
  use rg

  integer          :: N, nit
  double precision :: lambda

  integer :: it ! iteration of the algorithm it goes from 1 to nit

  ! DEBUGGING PURPOSES
  logical :: DEBUG
  integer :: iostat

  ! of each iteration...
  double precision, dimension(:,:), allocatable :: HN     ! first hamiltonian
  double precision, dimension(:,:), allocatable :: H2N    ! doubled hamiltonian
  double precision, dimension(:,:), allocatable :: P      !

  double precision, dimension(:), allocatable :: evls
  double precision                            :: oldevl
  integer*8                                   :: sizeofspace

  character(20) :: folder, file

  DEBUG = .TRUE.
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
  print*, "+------------------------------------------+"

  write(*,"(A)",advance='no') " + Type N, nit, lambda, folder, file: "
  read (*,*, iostat=iostat) N, nit, lambda, folder, file

  call system("mkdir -p ./data/"//trim(folder))
  open(1, file = "./data/"//trim(folder)//"/"//trim(file)//".csv")

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

  print*, "+ Iterating...                             +"
  ! iterations
  do it = 1, nit, 1
    if(modulo(it*10,nit)==0) then
      write(*,'(A,I6,A,I6,A)',advance='yes') " + ", it, "  /", nit, "                          +"
    end if

    allocate(H2N(2**(2*N),2**(2*N)))
    sizeofspace = 2*sizeofspace

    if(it==1) then
      H2N = mat_tensor_I(HN) + I_tensor_mat(HN) + init_interaction_H(N)
    else
      H2N = mat_tensor_I(HN) + I_tensor_mat(HN) + matmul(matmul(transpose(P),init_interaction_H(N)),P)
      deallocate(P)
    end if

    call findkEigenvalue(H2N, evls, 2**N, P)

    HN = matmul(matmul(transpose(P),H2N),P)
 
    if(abs(oldevl - evls(1)/sizeofspace)<1d-10) then
      print*, "+ Algorithm has reached convergence"
      !stop
    end if
    
    write(1,*) evls(1)/sizeofspace
    oldevl = evls(1)/sizeofspace
    
    

    deallocate(H2N,evls)

  end do
  print*, "+------------------------------------------+"
  print*, "+ Done!                                    +"
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program ising_rg
